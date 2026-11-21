def variant_key(rec):
    alts = tuple(rec.alts) if rec.alts else ()
    return (rec.chrom, rec.pos, rec.ref, alts)


def chrom_sort_key(chrom):
    c = chrom.replace("chr", "")
    order = {str(i): i for i in range(1, 23)}
    order.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    return order.get(c, 1000), chrom


def get_filter_string(rec):
    filt = list(rec.filter.keys())
    return "." if not filt else ";".join(filt)


def get_qual_string(rec):
    return "." if rec.qual is None else str(rec.qual)


def infer_tumor_normal_names(vcf):
    sample_names = list(vcf.header.samples)
    tumor = tumor_sample_name if tumor_sample_name in sample_names else None
    normal = normal_sample_name if normal_sample_name in sample_names else None

    if tumor is None:
        for s in sample_names:
            if "tumor" in s.lower() or "tumour" in s.lower():
                tumor = s
                break
    if normal is None:
        for s in sample_names:
            if "normal" in s.lower() or "control" in s.lower():
                normal = s
                break

    if tumor is None and len(sample_names) >= 1:
        tumor = sample_names[0]
    if normal is None and len(sample_names) >= 2:
        normal = sample_names[1]

    return tumor, normal


def _as_int_list(value):
    if value is None:
        return None
    if isinstance(value, str):
        if value in (".", ""):
            return None
        parts = value.split(",")
        out = []
        for p in parts:
            if p in (".", ""):
                return None
            out.append(int(float(p)))
        return out
    if isinstance(value, (tuple, list)):
        out = []
        for x in value:
            if x is None or x == ".":
                return None
            out.append(int(x))
        return out
    return [int(value)]


def _get_sample_field(sample, key):
    try:
        return sample.get(key)
    except Exception:
        return None


def extract_ad_dp_from_sample(rec, sample_name, caller):
    """
    Return:
        ad_tuple: (ref_count, alt_count1, ...)
        dp_int
    """
    if sample_name is None or sample_name not in rec.samples:
        return None, None

    sample = rec.samples[sample_name]
    ref = rec.ref
    alts = list(rec.alts or [])
    n_alt = len(alts)

    # 1) 通用 AD + DP
    ad = _as_int_list(_get_sample_field(sample, "AD"))
    dp = _get_sample_field(sample, "DP")
    dp = None if dp in (None, ".") else int(dp)
    if ad is not None and len(ad) >= 1 + n_alt:
        return tuple(ad[:1 + n_alt]), dp if dp is not None else int(sum(ad))

    caller_l = caller.lower()

    # 2) VarScan2: RD + AD
    if caller_l.startswith("varscan"):
        rd = _get_sample_field(sample, "RD")
        ad_alt = _get_sample_field(sample, "AD")
        try:
            if rd not in (None, ".") and ad_alt not in (None, "."):
                ad2 = (int(rd), int(ad_alt))
                return ad2, int(sum(ad2))
        except Exception:
            pass

    # 3) SomaticSniper: DP4 = refF,refR,altF,altR
    dp4 = _get_sample_field(sample, "DP4")
    dp4 = _as_int_list(dp4)
    if dp4 is not None and len(dp4) >= 4:
        ref_count = dp4[0] + dp4[1]
        alt_count = dp4[2] + dp4[3]
        ad2 = (ref_count, alt_count)
        return ad2, int(sum(ad2))

    # 4) Strelka2 indel: TAR/TIR
    tar = _as_int_list(_get_sample_field(sample, "TAR"))
    tir = _as_int_list(_get_sample_field(sample, "TIR"))
    if tar is not None and tir is not None and n_alt == 1:
        ad2 = (tar[0], tir[0])
        return ad2, int(sum(ad2))

    # 5) Strelka2 SNV: refU / altU
    if n_alt == 1:
        ref_tag = f"{ref}U"
        alt_tag = f"{alts[0]}U"
        refu = _as_int_list(_get_sample_field(sample, ref_tag))
        altu = _as_int_list(_get_sample_field(sample, alt_tag))
        if refu is not None and altu is not None:
            ad2 = (refu[0], altu[0])
            return ad2, int(sum(ad2))

    # 6) 只有 DP
    if dp is not None:
        return None, dp

    return None, None


def pad_ad(ad_tuple, target_len):
    if ad_tuple is None:
        return None
    if len(ad_tuple) >= target_len:
        return tuple(ad_tuple[:target_len])
    return tuple(list(ad_tuple) + [0] * (target_len - len(ad_tuple)))


def max_ad(ad_list, allele_count):
    vals = [pad_ad(ad, allele_count) for ad in ad_list if ad is not None]
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    out = [0] * allele_count
    for i in range(allele_count):
        out[i] = max(v[i] for v in vals)
    return tuple(out)


def max_dp(dp_list):
    vals = [int(v) for v in dp_list if v not in (None, ".")]
    return max(vals) if vals else None


def pick_gt(sample_payloads):
    for payload in sample_payloads:
        gt = payload.get("GT")
        if gt is not None:
            return gt
    return (None, None)


def create_output_record(out_vcf, src_rec):
    filt = list(src_rec.filter.keys())
    if not filt:
        filt = None

    new_rec = out_vcf.new_record(
        contig=src_rec.chrom,
        start=src_rec.start,
        stop=src_rec.stop,
        id=src_rec.id,
        alleles=(src_rec.ref,) + tuple(src_rec.alts or ()),
        qual=src_rec.qual,
    )

    # for f in src_rec.filter.keys():
    #     new_rec.filter.add(f)

    # 尽量保留原 INFO（前提是 header 里存在）
    for k in src_rec.info.keys():
        if k in out_vcf.header.info and k not in {
            "CALLER", "CALLERS", "NCALLERS",
            "CALLER_FILTERS", "CALLER_QUALS",
            "CONSENSUS_PASS", "TUMOR_SAMPLE", "NORMAL_SAMPLE"
        }:
            try:
                new_rec.info[k] = src_rec.info[k]
            except Exception:
                pass

    return new_rec
