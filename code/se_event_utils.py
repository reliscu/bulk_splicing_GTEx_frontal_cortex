from collections import defaultdict

def find_overlapping_cds(exon_start, exon_end, transcript_cds_rows):
    """Return the CDS row that overlaps this exon, or None."""
    for _, cds in transcript_cds_rows.iterrows():
        if cds['end'] >= exon_start and cds['start'] <= exon_end:
            return cds
    return None

def classify_overlap(exon_start, exon_end, cds_row):
    if cds_row['start'] == exon_start and cds_row['end'] == exon_end:
        return "fully_coding"
    elif cds_row['start'] >= exon_start and cds_row['end'] <= exon_end:
        return "partially_coding" # CDS subset of exon (should only happen when part of exon is UTR)
    else:
        return "unexpected"    # one side hangs over — rare
    
def _flanks_match(exons, es, ee, us_intron_start, ds_intron_end):
    left = [x for x in exons if x['end'] < es]
    right = [x for x in exons if x['start'] > ee]
    if not left or not right:
        return False                                  # terminal exon in this transcript
    left_exon = max(left, key=lambda x: x['end'])
    right_exon = min(right, key=lambda x: x['start'])
    return (left_exon['end'] + 1 == us_intron_start) and (right_exon['start'] - 1 == ds_intron_end)

def classify_isoform_detailed(exons, es, ee, us_intron_start, ds_intron_end, strand):
    """Returns (bucket, detail). detail carries the boundary sub-type for
    exon_diff_boundary, else None."""
    exact = overlap = None
    has_up = has_down = exonic_overlap = False
    for ex in exons:
        s, e = ex['start'], ex['end']
        if e < es:
            has_up = True
        if s > ee:
            has_down = True
        if s <= ee and e >= es:
            exonic_overlap = True
            if s == es and e == ee:
                exact = ex
            else:
                overlap = ex
    if exact is not None:
        s, e = exact['start'], exact['end']   # use exact, not last loop var
        bucket = ("compatible" if _flanks_match(exons, es, ee, us_intron_start, ds_intron_end)
                else "exon_diff_junction")
        return bucket, {'start': s, 'end': e}
    if overlap is not None:
        s, e = overlap['start'], overlap['end']
        if s == es:                                   # shares genomic-left boundary
            kind = "alt_5ss" if strand == '+' else "alt_3ss"
        elif e == ee:                                 # shares genomic-right boundary
            kind = "alt_3ss" if strand == '+' else "alt_5ss"
        else:
            kind = "overlapping_exon"
        return "exon_diff_boundary", {'start': s, 'end': e, 'kind': kind}
    if has_up and has_down and not exonic_overlap:
        return "exon_skipped", None
    return "locus_not_covered", None

def map_exon_to_protein(es, ee, cds_obj, strict=True):
    """Returns the CDS/aa annotation dict, or None if the transcript is
    noncoding or the exon falls in UTR. CDS rows MUST be in translation order."""
    cds = _cds_rows(cds_obj)
    if cds is None:
        return None                                   # noncoding transcript
    overlapping = next((c for c in cds if c['end'] >= es and c['start'] <= ee), None)
    if overlapping is None:
        return None                                   # exon is UTR here
    exon_cds_start, exon_cds_end = overlapping['start'], overlapping['end']
 
    cds_offset, gtf_frame = 0, None
    for c in cds:                                     # translation order
        if c['start'] == exon_cds_start and c['end'] == exon_cds_end:
            gtf_frame = c['frame']
            break
        cds_offset += c['end'] - c['start'] + 1
    rel_exon_start = cds_offset
    rel_exon_end = cds_offset + (exon_cds_end - exon_cds_start)
 
    first_frame = cds[0]['frame']                     # row 0 == translation-start CDS
    this_cds_frame = (first_frame - rel_exon_start) % 3
    if this_cds_frame != gtf_frame:
        msg = f"frame mismatch: computed {this_cds_frame}, GTF {gtf_frame} (check CDS ordering)"
        if strict:
            raise AssertionError(msg)
        return {'error': msg}
 
    coding_nt_length = exon_cds_end - exon_cds_start + 1
    overlap_type = ("fully_coding" if exon_cds_start == es and exon_cds_end == ee
                    else "partially_coding" if exon_cds_start >= es and exon_cds_end <= ee
                    else "unexpected")
    return {
        'aa_start': rel_exon_start // 3,
        'aa_end': rel_exon_end // 3,
        'exon_cds_start': exon_cds_start,
        'exon_cds_end': exon_cds_end,
        'coding_nt_length': coding_nt_length,
        'overlap_type': overlap_type,
        'gtf_frame': gtf_frame,
        'clean_start': gtf_frame == 0,
        'clean_end': (coding_nt_length - gtf_frame) % 3 == 0,
        'frame_preserving': coding_nt_length % 3 == 0,
        
    }
 
def _exon_rows(obj):
    if obj is None:
        return []
    if hasattr(obj, "iterrows"):
        return [{'start': int(r['start']), 'end': int(r['end'])} for _, r in obj.iterrows()]
    return [{'start': int(e['start']), 'end': int(e['end'])} for e in obj]

def annotate_event(event, transcripts_by_gene, exons_by_transcript, cds_by_transcript,
                   strict=True):
    """event: {'chr','strand','es','ee','gene','us_intron_start','ds_intron_end'}."""
    rec = {'meta': dict(event), 'cluster_id': None,
           'compatible': {}, 'exon_diff_junction': {},
           'exon_diff_boundary': {}, 'exon_skipped': {}, 'locus_not_covered': {}}
    es, ee = event['es'], event['ee']
    us, ds = event['us_intron_start'], event['ds_intron_end']
    
    for _, row in transcripts_by_gene.get(event['gene'], []).iterrows():
        t = row['transcript']
        ttype = row.get('transcript_type') 
        ttag = row.get('tag')
        exons = _exon_rows(exons_by_transcript.get(t))
        if not exons:
            continue
        
        bucket, detail = classify_isoform_detailed(exons, es, ee, us, ds, event['strand']) 
        if bucket in ("compatible", "exon_diff_junction"):
            mapped = map_exon_to_protein(es, ee, cds_by_transcript.get(t), strict=strict)
            entry = mapped if mapped is not None else {'overlap_type': 'noncoding_or_utr'}
        elif bucket == "exon_diff_boundary":
            entry = dict(detail)
            sib_es, sib_ee = detail['start'], detail['end']
            mapped = map_exon_to_protein(sib_es, sib_ee, cds_by_transcript.get(t), strict=strict)
            if mapped is not None:
                entry.update(mapped)
            else:
                entry['overlap_type'] = 'noncoding_or_utr'
        else:
            entry = {}
            
        entry['transcript_type'] = ttype
        entry['transcript_tag'] = ttag
        rec[bucket][t] = entry
        
    return rec

def _cds_rows(obj):
    if obj is None:
        return None
    if hasattr(obj, "iterrows"):
        return [{'start': int(r['start']), 'end': int(r['end']), 'frame': int(r['frame'])}
                for _, r in obj.iterrows()]
    return [{'start': int(c['start']), 'end': int(c['end']), 'frame': int(c['frame'])} for c in obj]

def cluster_events(events):
    """events: list of dicts with chr,strand,es,ee; adds 'cluster_id' in place."""
    by_csg = defaultdict(list)
    for ev in events:
        by_csg[(ev['chr'], ev['strand'], ev['gene'])].append(ev)
    cid = 0
    for grp in by_csg.values():
        grp.sort(key=lambda x: (x['es'], x['ee']))
        cur_end = None
        for ev in grp:
            if cur_end is None or ev['es'] > cur_end:
                cid += 1
                cur_end = ev['ee']
            else:
                cur_end = max(cur_end, ev['ee'])
            ev['cluster_id'] = cid
    return events

def mark_sibling_variants(event_info, event_coords):
    """Flag whether a SE variant matches a *called event*"""
    coords_by_cluster = defaultdict(set)
    for ev, cid in event_coords.items():
        coords_by_cluster[cid[0]].add((cid[1], cid[2]))   # cid = (cluster_id, es, ee)
    
    # all transcripts that appear as compatible in any event
    compatible_transcripts = set()
    for rec in event_info.values():
        compatible_transcripts.update(rec['compatible'].keys())
        
    for rec in event_info.values():
        cid = rec['cluster_id']
        cluster = coords_by_cluster.get(cid, set())
        parent_es = rec['meta']['es']
        parent_ee = rec['meta']['ee']
        
        # boundary variants: same cluster + shares at least one exon boundary
        for t, sib in rec['exon_diff_boundary'].items():
            sib_es, sib_ee = sib['start'], sib['end']
            sib['is_called_sibling'] = (
                (sib_es, sib_ee) in cluster
                and (sib_es == parent_es or sib_ee == parent_ee)
            )
            
        # junction variants: same exon coords, different flanking partner
        for t, d in rec['exon_diff_junction'].items():
            if not isinstance(d, dict):
                continue
            d['is_called_sibling'] = (parent_es, parent_ee) in cluster
            
        # skipped variants
        for t, d in rec['exon_skipped'].items():
            if not isinstance(d, dict):
                continue
            d['is_called_sibling'] = t in compatible_transcripts
            
    return event_info

def score_transcript(t, x):
    tag = str(x['transcript_tag'])
    return (
        int("MANE_Select" in tag),
        int("MANE_Plus_Clinical" in tag),
        int("appris_principal" in tag),
        int("basic" in tag),
        int("CCDS" in tag),
        int("GENCODE_Primary" in tag),
        x['coding_nt_length'] if "coding_nt_length" in x.keys() else 0, # prefer longer (more complete) CDS
        t # deterministic alphabetical tiebreak
    ) 

def _cds_rows(obj):
    if obj is None:
        return None
    if hasattr(obj, "iterrows"):
        return [{'start': int(r['start']), 'end': int(r['end']), 'frame': int(r['frame'])}
                for _, r in obj.iterrows()]
    return [{'start': int(c['start']), 'end': int(c['end']), 'frame': int(c['frame'])} for c in obj]
 
def get_coding_nt_length(transcript, cds_by_transcript):
    cds = cds_by_transcript.get(transcript)
    if cds is None:
        return 0
    return sum(c['end'] - c['start'] + 1 for c in _cds_rows(cds))