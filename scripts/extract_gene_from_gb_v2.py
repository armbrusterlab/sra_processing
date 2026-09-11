from pathlib import Path
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


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


def extract_gene_to_gb(gb_file, outname, overwrite=True, asm=None, target=None, 
                       target_type=None, buffer_upstream=0, buffer_downstream=0, 
                       translate=False, mask_file=None, include_pseudogenes=True):
    '''
    Version that writes to GenBank with CDS annotations.
    
    Parameters:
    -----------
    gb_file : str
        Path to input GenBank file
    outname : str
        Path to output GenBank file with extracted genes
    overwrite : bool
        Whether to overwrite existing output files
    asm : str
        Assembly accession (if None, uses parent directory name)
    target : list
        List of target values to filter by
    target_type : str
        Qualifier type to filter on (e.g., "gene", "product")
    buffer_upstream : int
        Number of upstream bases to include as buffer
    buffer_downstream : int
        Number of downstream bases to include as buffer
    translate : bool
        Whether to translate CDS sequences to protein
    mask_file : str
        Optional path to output masked GenBank file
    include_pseudogenes : bool
        Whether to include pseudogenes (with /pseudo tag) in output
    '''
    assembly_accession = asm if asm is not None else Path(gb_file).parent.name
    
    if overwrite:
        with open(outname, "w"):
            pass
        if mask_file:
            with open(mask_file, "w"):
                pass

    for record in SeqIO.parse(gb_file, "genbank"):
        # For mask_file: create a mutable copy of the sequence
        if mask_file:
            masked_seq = list(str(record.seq))
            regions_to_mask = []
        
        for feature in record.features:
            # Check if this is a CDS or pseudogene
            is_cds = (feature.type == "CDS")
            is_pseudo = feature.qualifiers.get("pseudo", [False])[0] if "pseudo" in feature.qualifiers else False
            
            # If include_pseudogenes is False, skip pseudogenes
            if not include_pseudogenes and is_pseudo:
                continue
            
            if is_cds or (include_pseudogenes and is_pseudo and feature.type == "CDS"):
                # For pseudogenes, we treat them similarly to CDS
                # but mark them differently in output
                
                if target is None or target_type is None or feature.qualifiers.get(target_type, ["N/A"])[0] in target:
                    loc = feature.location
                    
                    # Store the region to mask (with buffers)
                    if mask_file:
                        # Get the full buffered region coordinates
                        if isinstance(loc, CompoundLocation):
                            for i, part in enumerate(loc.parts):
                                is_first = (i == 0)
                                is_last = (i == len(loc.parts) - 1)
                                
                                start = part.start - (buffer_upstream if is_first else 0)
                                end = part.end + (buffer_downstream if is_last else 0)
                                start = max(0, start)
                                end = min(len(record.seq), end)
                                regions_to_mask.append((start, end))
                        else:
                            start = max(0, loc.start - buffer_upstream)
                            end = min(len(record.seq), loc.end + buffer_downstream)
                            regions_to_mask.append((start, end))
                    
                    # Extract the buffered sequence (handles join/complement)
                    flanked_seq, cds_start_in_buffer, cds_end_in_buffer = extract_with_buffer(
                        record, loc, buffer_upstream, buffer_downstream
                    )

                    # important: strand info is encoded at feature.location.strand. 
                    # If negative strand, need to reverse complement.
                    if loc.strand == -1:
                        flanked_seq = flanked_seq.reverse_complement()

                    if translate and not is_pseudo:
                        # Only translate true CDS, not pseudogenes
                        flanked_seq = flanked_seq.translate(to_stop=True)
                    elif translate and is_pseudo:
                        # For pseudogenes, don't translate but note it
                        pass

                    locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                    product = feature.qualifiers.get("product", ["N/A"])[0]
                    gene = feature.qualifiers.get("gene", ["N/A"])[0]
                    
                    # Get pseudogene-specific annotations if available
                    pseudo_note = feature.qualifiers.get("note", ["N/A"])[0] if is_pseudo else None

                    # Compute CDS coordinates within the buffered sequence
                    cds_start_in_buffer = buffer_upstream
                    cds_end_in_buffer = buffer_upstream + len(feature.location)

                    # Create SeqRecord for GenBank output
                    buffer_string = f" | Sequence with upstream buffer <= {buffer_upstream} bp, downstream buffer <= {buffer_downstream} bp" if (buffer_upstream != 0 or buffer_downstream != 0) else ""
                    pseudo_tag = " [PSEUDOGENE]" if is_pseudo else ""
                    new_record = SeqRecord(
                        flanked_seq,
                        id=f"{assembly_accession}-{protein_id}-{locus_tag}",
                        name=locus_tag,
                        description=f"{gene} | {product}{pseudo_tag}{buffer_string}"
                    )

                    new_record.annotations["molecule_type"] = "protein" if translate and not is_pseudo else "DNA"
                    if is_pseudo:
                        new_record.annotations["pseudo"] = True
                    new_record.seq = flanked_seq

                    # Build CDS feature (or pseudo-CDS feature)
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
                    
                    # Add pseudo qualifier if it's a pseudogene
                    if is_pseudo:
                        cds_feature.qualifiers["pseudo"] = [""]  # /pseudo tag with no value
                        if pseudo_note and pseudo_note != "N/A":
                            cds_feature.qualifiers["note"] = [pseudo_note]

                    new_record.features.append(cds_feature)

                    with open(outname, "a") as output_handle:
                        SeqIO.write(new_record, output_handle, "genbank")

        # After processing all features, write the masked GenBank file for this record
        if mask_file and regions_to_mask:
            # Merge overlapping regions
            regions_to_mask.sort()
            merged_regions = []
            for start, end in regions_to_mask:
                if not merged_regions or start > merged_regions[-1][1]:
                    merged_regions.append([start, end])
                else:
                    merged_regions[-1][1] = max(merged_regions[-1][1], end)
            
            # Apply masking
            for start, end in merged_regions:
                for pos in range(start, min(end, len(masked_seq))):
                    masked_seq[pos] = 'N'
            
            # Create masked record
            masked_record = SeqRecord(
                Seq(''.join(masked_seq)),
                id=record.id,
                name=record.name,
                description=record.description,
                annotations=record.annotations.copy()
            )
            
            # Copy only features that don't overlap with masked regions
            for feature in record.features:
                keep_feature = True
                if hasattr(feature.location, 'start') and hasattr(feature.location, 'end'):
                    feat_start = feature.location.start
                    feat_end = feature.location.end
                    # Check if feature overlaps with any masked region
                    for start, end in merged_regions:
                        if not (feat_end <= start or feat_start >= end):
                            keep_feature = False
                            break
                
                if keep_feature:
                    masked_record.features.append(feature)
            
            # Write masked record
            with open(mask_file, "a") as mask_handle:
                SeqIO.write(masked_record, mask_handle, "genbank")


# Example usage:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to extract selected genes with buffer and save them as GB format for the breseq step of the pipeline.")

    # required
    parser.add_argument("-i", "--gb_file", type=str, help='Path to GenBank file from which selected genes should be extracted.')
    parser.add_argument("-o", "--outname", type=str, help='Name of output .gb format file.')
    parser.add_argument("-g", "--target", nargs='*', help='List of genes to extract.')
    parser.add_argument("-t", "--target_type", type=str, help='Type of target (e.g. "locus_tag", "gene", or "product").')

    # optional
    parser.add_argument("-a", "--asm", type=str, default=None, help='Assembly accession.')
    parser.add_argument("-u", "--buffer_upstream", type=int, default=900, help='Number of upstream bases to include as buffer around each gene.')
    parser.add_argument("-d", "--buffer_downstream", type=int, default=900, help='Number of downstream bases to include as buffer around each gene.')

    args = parser.parse_args()
    extract_gene_to_gb(gb_file=args.gb_file, outname=args.outname, asm=args.asm, target=args.target, target_type=args.target_type,
                       buffer_upstream=args.buffer_upstream, buffer_downstream=args.buffer_downstream)