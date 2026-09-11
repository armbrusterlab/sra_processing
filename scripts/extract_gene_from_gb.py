from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# note: extract_gene doesn't yet account for compound locations like extract_gene_to_gb does
def extract_gene(gb_file, outname, overwrite = True, asm = None, target = None, target_type = None, buffer_upstream=0, buffer_downstream=0, translate=False):
    '''
    gb_file = path to GenBank.
    To outname, write a FASTA either containing just the gene, or if either target or target_type are None, all genes in FASTA.
    If overwrite, clears the contents of outname; otherwise continue to build upon FASTA located at outname.

    asm: assembly accession (i.e. name of dir to look in when searching genome_db for GBFF to retrieve metadata from later); a failsafe for inference from gb_file path

    target = query name (list); target_type indicates where in the record to search for the target
    (e.g. target = ["mucA", "wspF"], target_type = "gene"; 
    target = ["NP_249454.1"], target_type = "protein_id"; 
    target = ["PA0763"], target_type = "locus_tag")
    Please note that if the target matches multiple sequences, then multiple sequences will be written to outname.

    Buffers are defined in terms of bp.

    If using translate==True, then please use buffer_upstream==0, as translation starts at the start location of the sequence.
    '''
    assembly_accession = asm if asm is not None else Path(gb_file).parent.name # even if it's not actually an assembly accession, want the dir name so we can use the fasta's seq id to find the original gb_file later
    # in the context of Nextflow this might be unreliable, but it can serve as a backup
    # it turns out that it's unnecessarily complicated to get the assembly accession from a GenBank file

    if overwrite:
        with open(outname, "w"):
            pass

    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if target is None or target_type is None or feature.qualifiers.get(target_type, ["N/A"])[0] in target: 
                    # expect the 0-th element to be the only element
                    # Get location, handling strand direction
                    loc = feature.location
                    # Extract start/end with flanking buffer
                    start = max(0, loc.start - buffer_upstream)
                    end = min(len(record.seq), loc.end + buffer_downstream)
                    # Slice and extract
                    flanked_seq = record.seq[start:end]

                    # important: strand info is encoded at feature.location.strand. If negative strand, need to reverse complement.
                    if loc.strand == -1:
                        flanked_seq = flanked_seq.reverse_complement()

                    if translate:
                        flanked_seq = flanked_seq.translate(to_stop=True)

                    locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                    product = feature.qualifiers.get("product", ["N/A"])[0]
                    gene = feature.qualifiers.get("gene", ["N/A"])[0]

                    # Create a new SeqRecord
                    new_record = SeqRecord(
                        flanked_seq,
                        id=f"{assembly_accession}-{protein_id}-{locus_tag} | {gene} | {product} |",
                        description=f"buffer: upstream +{loc.start-start} bp, downstream +{end-loc.end} bp"
                    )

                    # Save to FASTA
                    with open(outname, "a") as output_handle:
                        SeqIO.write(new_record, output_handle, "fasta")


def extract_gene_to_gb(gb_file, outname, overwrite = True, asm = None, target = None, target_type = None, buffer_upstream=0, buffer_downstream=0, translate=False):
    '''
    Version that writes to GenBank with CDS annotations.
    '''
    assembly_accession = asm if asm is not None else Path(gb_file).parent.name # even if it's not actually an assembly accession, want the dir name so we can use the fasta's seq id to find the original gb_file later
    # in the context of Nextflow this might be unreliable, but it can serve as a backup
    # it turns out that it's unnecessarily complicated to get the assembly accession from a GenBank file

    if overwrite:
        with open(outname, "w"):
            pass

    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if target is None or target_type is None or feature.qualifiers.get(target_type, ["N/A"])[0] in target: 
                    # expect the 0-th element to be the only element
                    # Get location, handling strand direction
                    loc = feature.location
                    # Extract buffered sequence (handles join/complement)
                    flanked_seq, cds_start_in_buffer, cds_end_in_buffer = extract_with_buffer(
                        record, loc, buffer_upstream, buffer_downstream
                    )

                    # important: strand info is encoded at feature.location.strand. If negative strand, need to reverse complement.
                    if loc.strand == -1:
                        flanked_seq = flanked_seq.reverse_complement()

                    if translate:
                        flanked_seq = flanked_seq.translate(to_stop=True)

                    locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                    product = feature.qualifiers.get("product", ["N/A"])[0]
                    gene = feature.qualifiers.get("gene", ["N/A"])[0]

                    # Compute CDS coordinates within the buffered sequence
                    cds_start_in_buffer = buffer_upstream
                    cds_end_in_buffer = buffer_upstream + len(feature.location)

                    # Create SeqRecord for GenBank output
                    buffer_string = f" | Sequence with upstream buffer <= {buffer_upstream} bp, downstream buffer <= {buffer_downstream} bp" if (buffer_upstream != 0 or buffer_downstream != 0) else ""
                    new_record = SeqRecord(
                        flanked_seq,
                        id=f"{assembly_accession}-{protein_id}-{locus_tag}",
                        name=locus_tag,
                        description=f"{gene} | {product}{buffer_string}" # finding exact buffer lengths is more complicated when considering compound locations
                    )

                    new_record.annotations["molecule_type"] = "protein" if translate else "DNA"

                    # Use the case set in extract_with_buffer
                    new_record.seq = flanked_seq

                    # Build CDS feature
                    # new_strand = loc.strand
                    new_strand = 1 # since the sequence has been reverse-complemented if it was on the negative strand, always write CDS as if on positive strand
                    cds_feature = SeqFeature(
                        FeatureLocation(cds_start_in_buffer, cds_end_in_buffer, strand=1),
                        type="CDS",
                        qualifiers={
                            "locus_tag": [locus_tag],
                            "gene": [gene],
                            "product": [product],
                            "protein_id": [protein_id]
                        }
                    )

                    # Attach feature
                    new_record.features.append(cds_feature)

                    # Write GenBank instead of FASTA
                    with open(outname, "a") as output_handle:
                        SeqIO.write(new_record, output_handle, "genbank") # doesn't write with case unfortunately
                        # write_genbank_with_case(new_record, output_handle) # should work, but comparatively risky, and the only benefit is that the GB is easier for the user to read when the GB is primarily meant as an input to breseq


def extract_with_buffer(record, loc, buffer_upstream, buffer_downstream):
    """
    Extract a buffered sequence for FeatureLocation or CompoundLocation.
    For joined features:
        - upstream buffer applies only to the first part
        - downstream buffer applies only to the last part
        - For negative strand, fragments are joined in reverse order
        - Buffer regions are lowercase, CDS regions are uppercase
    """
    seq = record.seq

    # Simple case: single interval
    if not isinstance(loc, CompoundLocation):
        # Extract the full buffered sequence (all lowercase initially)
        start = max(0, loc.start - buffer_upstream)
        end = min(len(seq), loc.end + buffer_downstream)
        full_seq = seq[start:end]
        
        # Convert to mutable list for case manipulation
        seq_list = list(str(full_seq))
        
        # Determine CDS region within the buffer
        cds_start_in_buffer = buffer_upstream
        cds_end_in_buffer = buffer_upstream + len(loc)
        
        # Uppercase the CDS region
        for pos in range(cds_start_in_buffer, min(cds_end_in_buffer, len(seq_list))):
            seq_list[pos] = seq_list[pos].upper()
        
        # Lowercase the buffer regions (already lowercase, but ensure it)
        for pos in range(0, cds_start_in_buffer):
            seq_list[pos] = seq_list[pos].lower()
        for pos in range(cds_end_in_buffer, len(seq_list)):
            seq_list[pos] = seq_list[pos].lower()
        
        buffered = Seq(''.join(seq_list))
        return buffered, cds_start_in_buffer, cds_end_in_buffer

    # CompoundLocation case
    parts = []
    part_boundaries = []  # Store (start_in_part, end_in_part, is_cds) for each part
    total_len = 0

    # Determine if we need to process fragments in reverse order
    strand = loc.strand
    part_indices = range(len(loc.parts)) if strand != -1 else range(len(loc.parts) - 1, -1, -1)

    for i in part_indices:
        part = loc.parts[i]
        is_first_in_original = (i == 0)
        is_last_in_original = (i == len(loc.parts) - 1)
        
        # Determine buffer application
        if strand == -1:
            is_first_in_processed = (i == len(loc.parts) - 1)
            is_last_in_processed = (i == 0)
            apply_upstream = is_first_in_processed
            apply_downstream = is_last_in_processed
        else:
            apply_upstream = is_first_in_original
            apply_downstream = is_last_in_original

        # Get the extracted part (all lowercase initially)
        start = part.start - (buffer_upstream if apply_upstream else 0)
        end = part.end + (buffer_downstream if apply_downstream else 0)
        start = max(0, start)
        end = min(len(seq), end)
        
        part_seq = seq[start:end]
        
        # Determine CDS boundaries within this part
        if apply_upstream:
            cds_start_in_part = buffer_upstream
        else:
            cds_start_in_part = 0
        
        if apply_downstream:
            cds_end_in_part = len(part_seq) - buffer_downstream
        else:
            cds_end_in_part = len(part_seq)
        
        # Convert to list for case manipulation
        part_list = list(str(part_seq))
        
        # Uppercase the CDS region within this part
        for pos in range(cds_start_in_part, min(cds_end_in_part, len(part_list))):
            part_list[pos] = part_list[pos].upper()
        
        # Lowercase buffer regions (already lowercase, but ensure it)
        for pos in range(0, cds_start_in_part):
            part_list[pos] = part_list[pos].lower()
        for pos in range(cds_end_in_part, len(part_list)):
            part_list[pos] = part_list[pos].lower()
        
        parts.append(Seq(''.join(part_list)))
        total_len += len(part)

    # Join all extracted parts
    buffered = parts[0].__class__("").join(parts)

    # CDS starts after upstream buffer on the first part
    cds_start = buffer_upstream
    cds_end = buffer_upstream + len(loc)

    return buffered, cds_start, cds_end

def write_genbank_with_case(record, handle):
    """
    Write GenBank file preserving sequence case.
    """
    # First, get the formatted GenBank string
    gb_str = record.format("genbank")
    lines = gb_str.split('\n')
    
    new_lines = []
    in_sequence = False
    seq_printed = False
    
    for line in lines:
        if line.startswith('ORIGIN'):
            in_sequence = True
            new_lines.append(line)
            seq_printed = False
        elif in_sequence and line.strip() and not line.startswith('//'):
            # Skip the original sequence lines
            continue
        elif in_sequence and line.startswith('//'):
            # Write our custom sequence with preserved case
            seq_str = str(record.seq)
            # Format sequence with line breaks (60 characters per line)
            for i in range(0, len(seq_str), 60):
                chunk = seq_str[i:i+60]
                position = i + 1
                new_lines.append(f"{position:9} {chunk}")
            new_lines.append('//')
            in_sequence = False
        else:
            new_lines.append(line)
    
    handle.write('\n'.join(new_lines) + '\n')
