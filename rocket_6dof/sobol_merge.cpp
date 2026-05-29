//  sobol_merge.cpp  --  Merge Sobol indices from multiple chunks
//
//  Given M chunk CSVs (each at base N=N_chunk, independent seeds),
//  produce a single merged CSV with point estimates from the chunk
//  average and uncertainty from chunk-to-chunk variance.
//
//  Per-chunk CSV columns (from sobol_runner.cpp):
//    output, input, S1, ST, S1_lo, S1_hi, ST_lo, ST_hi,
//    clean_share, clean_share_lo, clean_share_hi,
//    output_mean, output_variance
//
//  Merge math:
//    S1_merged   = (1/M) * sum_m S1_m
//    SE          = sqrt(var(S1_m across chunks) / M)
//    CI bounds   = S1_merged +/- 1.96 * SE  (normal approx)
//
//  This is statistically valid when chunks are independent samples
//  (different RNG seeds), which the chunk launch script guarantees.
//
//  Skips inputs where any chunk reports non-finite.  Writes one
//  output row per (output, input) tuple, same format as a single
//  Sobol run -- so the existing viewer can ingest it directly.
//
//  Usage:
//    sobol_merge OUT.csv  chunk0.csv chunk1.csv chunk2.csv chunk3.csv
//
//  Optional second pass: merge pair indices.  Use --pairs INPUTS_CSV
//  flag with chunks's pair CSVs:
//    sobol_merge --pairs PAIRS_OUT.csv  chunk0_pairs.csv ... chunk3_pairs.csv

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <string>
#include <map>

namespace {

// Strip leading/trailing whitespace
std::string strip(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a]==' '||s[a]=='\t'||s[a]=='\r'||s[a]=='\n')) ++a;
    while (b > a && (s[b-1]==' '||s[b-1]=='\t'||s[b-1]=='\r'||s[b-1]=='\n')) --b;
    return s.substr(a, b-a);
}

std::vector<std::string> split_comma(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') { out.push_back(strip(cur)); cur.clear(); }
        else { cur.push_back(c); }
    }
    out.push_back(strip(cur));
    return out;
}

// Parse a double; return false for "nan", "-nan", non-numeric
bool parse_double(const std::string& s, double& out) {
    if (s.empty()) return false;
    std::string lo;
    for (char c : s) lo.push_back(std::tolower(c));
    if (lo == "nan" || lo == "-nan" || lo == "inf" || lo == "-inf") return false;
    char* endp = nullptr;
    out = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return false;
    if (!std::isfinite(out)) return false;
    return true;
}

// Per-row data from a chunk CSV
struct RowKey {
    std::string output;
    std::string input;
    bool operator<(const RowKey& o) const {
        if (output != o.output) return output < o.output;
        return input < o.input;
    }
};
struct RowVals {
    std::vector<double> S1;          // one entry per chunk
    std::vector<double> ST;
    std::vector<double> clean_share;
    std::vector<double> output_mean;
    std::vector<double> output_var;
};

// Header column indices
struct ColIdx {
    int output    = -1;
    int input     = -1;
    int S1        = -1;
    int ST        = -1;
    int clean     = -1;
    int out_mean  = -1;
    int out_var   = -1;
};

ColIdx find_columns(const std::vector<std::string>& hdr) {
    ColIdx c;
    for (size_t i = 0; i < hdr.size(); ++i) {
        const auto& h = hdr[i];
        if (h == "output")               c.output   = (int)i;
        else if (h == "input")           c.input    = (int)i;
        else if (h == "S1")              c.S1       = (int)i;
        else if (h == "ST")              c.ST       = (int)i;
        else if (h == "clean_share")     c.clean    = (int)i;
        else if (h == "output_mean")     c.out_mean = (int)i;
        else if (h == "output_variance") c.out_var  = (int)i;
    }
    return c;
}

// Pair CSV columns: output, input_i, input_j, S_ij, S_ij_lo, S_ij_hi, output_mean, output_variance
struct PairKey {
    std::string output;
    std::string input_i;
    std::string input_j;
    bool operator<(const PairKey& o) const {
        if (output != o.output)   return output < o.output;
        if (input_i != o.input_i) return input_i < o.input_i;
        return input_j < o.input_j;
    }
};
struct PairVals {
    std::vector<double> Sij;
    std::vector<double> output_mean;
    std::vector<double> output_var;
};

// Mean of finite values, returns NaN if none
double mean_finite(const std::vector<double>& v) {
    double s = 0; int n = 0;
    for (double x : v) if (std::isfinite(x)) { s += x; ++n; }
    return n > 0 ? s/n : std::nan("");
}

// Sample SE = sqrt( var / n ).  Returns 0 if n<2.
double se_finite(const std::vector<double>& v) {
    int n = 0;
    double s = 0, s2 = 0;
    for (double x : v) if (std::isfinite(x)) { s += x; s2 += x*x; ++n; }
    if (n < 2) return 0.0;
    double mu = s/n;
    double var = s2/n - mu*mu;
    if (var < 0) var = 0;
    // Use n-1 unbiased variance estimator, then SE = sqrt(var_unb/n)
    double var_unb = var * (double)n / (double)(n-1);
    return std::sqrt(var_unb / (double)n);
}

int merge_first_order(const char* out_path,
                      const std::vector<std::string>& chunk_paths)
{
    std::map<RowKey, RowVals> rows;
    int n_chunks = 0;

    for (const auto& path : chunk_paths) {
        FILE* fp = std::fopen(path.c_str(), "r");
        if (!fp) { std::fprintf(stderr, "ERROR: cannot open %s\n", path.c_str()); return 1; }
        ++n_chunks;
        char buf[4096];
        if (!std::fgets(buf, sizeof(buf), fp)) { std::fclose(fp); continue; }
        std::vector<std::string> hdr = split_comma(buf);
        ColIdx col = find_columns(hdr);
        if (col.S1 < 0 || col.ST < 0 || col.output < 0 || col.input < 0) {
            std::fprintf(stderr, "ERROR: %s missing required columns\n", path.c_str());
            std::fclose(fp);
            return 1;
        }
        while (std::fgets(buf, sizeof(buf), fp)) {
            std::vector<std::string> r = split_comma(buf);
            if ((int)r.size() <= col.ST) continue;
            RowKey k{ r[col.output], r[col.input] };
            double v;
            auto& vals = rows[k];
            vals.S1.push_back(parse_double(r[col.S1], v) ? v : std::nan(""));
            vals.ST.push_back(parse_double(r[col.ST], v) ? v : std::nan(""));
            vals.clean_share.push_back(
                col.clean >= 0 && col.clean < (int)r.size()
                    && parse_double(r[col.clean], v) ? v : std::nan(""));
            vals.output_mean.push_back(
                col.out_mean >= 0 && col.out_mean < (int)r.size()
                    && parse_double(r[col.out_mean], v) ? v : std::nan(""));
            vals.output_var.push_back(
                col.out_var >= 0 && col.out_var < (int)r.size()
                    && parse_double(r[col.out_var], v) ? v : std::nan(""));
        }
        std::fclose(fp);
    }

    FILE* fo = std::fopen(out_path, "w");
    if (!fo) { std::fprintf(stderr, "ERROR: cannot write %s\n", out_path); return 1; }
    std::fprintf(fo,
        "output,input,S1,ST,S1_lo,S1_hi,ST_lo,ST_hi,"
        "clean_share,clean_share_lo,clean_share_hi,"
        "output_mean,output_variance,n_chunks,S1_se,ST_se\n");

    for (const auto& kv : rows) {
        const auto& k = kv.first;
        const auto& v = kv.second;
        double mS1 = mean_finite(v.S1);
        double mST = mean_finite(v.ST);
        double mCS = mean_finite(v.clean_share);
        double mOM = mean_finite(v.output_mean);
        double mOV = mean_finite(v.output_var);
        double seS1 = se_finite(v.S1);
        double seST = se_finite(v.ST);
        double seCS = se_finite(v.clean_share);
        // 95% CI = mean +/- 1.96 * SE
        double S1_lo = mS1 - 1.96*seS1, S1_hi = mS1 + 1.96*seS1;
        double ST_lo = mST - 1.96*seST, ST_hi = mST + 1.96*seST;
        double CS_lo = mCS - 1.96*seCS, CS_hi = mCS + 1.96*seCS;
        // Count valid chunks
        int n_valid = 0;
        for (double x : v.S1) if (std::isfinite(x)) ++n_valid;
        std::fprintf(fo, "%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%.6g,%.6g\n",
            k.output.c_str(), k.input.c_str(),
            mS1, mST, S1_lo, S1_hi, ST_lo, ST_hi,
            mCS, CS_lo, CS_hi, mOM, mOV,
            n_valid, seS1, seST);
    }
    std::fclose(fo);

    std::printf("merged %d chunks into %s (%zu rows)\n",
                n_chunks, out_path, rows.size());
    return 0;
}

int merge_pairs(const char* out_path,
                const std::vector<std::string>& chunk_paths)
{
    std::map<PairKey, PairVals> rows;
    int n_chunks = 0;
    for (const auto& path : chunk_paths) {
        FILE* fp = std::fopen(path.c_str(), "r");
        if (!fp) { std::fprintf(stderr, "ERROR: cannot open %s\n", path.c_str()); return 1; }
        ++n_chunks;
        char buf[4096];
        if (!std::fgets(buf, sizeof(buf), fp)) { std::fclose(fp); continue; }
        std::vector<std::string> hdr = split_comma(buf);
        int c_out = -1, c_in_i = -1, c_in_j = -1, c_Sij = -1, c_om = -1, c_ov = -1;
        for (size_t i = 0; i < hdr.size(); ++i) {
            if (hdr[i] == "output")               c_out  = (int)i;
            else if (hdr[i] == "input_i")         c_in_i = (int)i;
            else if (hdr[i] == "input_j")         c_in_j = (int)i;
            else if (hdr[i] == "S_ij")            c_Sij  = (int)i;
            else if (hdr[i] == "output_mean")     c_om   = (int)i;
            else if (hdr[i] == "output_variance") c_ov   = (int)i;
        }
        if (c_out < 0 || c_in_i < 0 || c_in_j < 0 || c_Sij < 0) {
            std::fprintf(stderr, "ERROR: %s missing pair columns\n", path.c_str());
            std::fclose(fp);
            return 1;
        }
        while (std::fgets(buf, sizeof(buf), fp)) {
            std::vector<std::string> r = split_comma(buf);
            if ((int)r.size() <= c_Sij) continue;
            PairKey k{ r[c_out], r[c_in_i], r[c_in_j] };
            double v;
            auto& vals = rows[k];
            vals.Sij.push_back(parse_double(r[c_Sij], v) ? v : std::nan(""));
            vals.output_mean.push_back(
                c_om >= 0 && c_om < (int)r.size() && parse_double(r[c_om], v) ? v : std::nan(""));
            vals.output_var.push_back(
                c_ov >= 0 && c_ov < (int)r.size() && parse_double(r[c_ov], v) ? v : std::nan(""));
        }
        std::fclose(fp);
    }

    FILE* fo = std::fopen(out_path, "w");
    if (!fo) { std::fprintf(stderr, "ERROR: cannot write %s\n", out_path); return 1; }
    std::fprintf(fo, "output,input_i,input_j,S_ij,S_ij_lo,S_ij_hi,output_mean,output_variance,n_chunks,S_ij_se\n");
    for (const auto& kv : rows) {
        const auto& k = kv.first;
        const auto& v = kv.second;
        double mSij = mean_finite(v.Sij);
        double seSij = se_finite(v.Sij);
        double mOM = mean_finite(v.output_mean);
        double mOV = mean_finite(v.output_var);
        double lo = mSij - 1.96*seSij;
        double hi = mSij + 1.96*seSij;
        int n_valid = 0;
        for (double x : v.Sij) if (std::isfinite(x)) ++n_valid;
        std::fprintf(fo, "%s,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%.6g\n",
            k.output.c_str(), k.input_i.c_str(), k.input_j.c_str(),
            mSij, lo, hi, mOM, mOV, n_valid, seSij);
    }
    std::fclose(fo);

    std::printf("merged %d chunks into %s (%zu rows)\n",
                n_chunks, out_path, rows.size());
    return 0;
}

// Group CSV columns: output, group, S1, ST, output_mean, output_variance
struct GroupKey {
    std::string output;
    std::string group;
    bool operator<(const GroupKey& o) const {
        if (output != o.output) return output < o.output;
        return group < o.group;
    }
};
struct GroupVals {
    std::vector<double> S1;
    std::vector<double> ST;
    std::vector<double> output_mean;
    std::vector<double> output_var;
};

int merge_groups(const char* out_path,
                 const std::vector<std::string>& chunk_paths)
{
    std::map<GroupKey, GroupVals> rows;
    int n_chunks = 0;
    for (const auto& path : chunk_paths) {
        FILE* fp = std::fopen(path.c_str(), "r");
        if (!fp) { std::fprintf(stderr, "ERROR: cannot open %s\n", path.c_str()); return 1; }
        ++n_chunks;
        char buf[4096];
        if (!std::fgets(buf, sizeof(buf), fp)) { std::fclose(fp); continue; }
        std::vector<std::string> hdr = split_comma(buf);
        int c_out = -1, c_g = -1, c_S1 = -1, c_ST = -1, c_om = -1, c_ov = -1;
        for (size_t i = 0; i < hdr.size(); ++i) {
            if (hdr[i] == "output")               c_out = (int)i;
            else if (hdr[i] == "group")           c_g   = (int)i;
            else if (hdr[i] == "S1")              c_S1  = (int)i;
            else if (hdr[i] == "ST")              c_ST  = (int)i;
            else if (hdr[i] == "output_mean")     c_om  = (int)i;
            else if (hdr[i] == "output_variance") c_ov  = (int)i;
        }
        if (c_out < 0 || c_g < 0 || c_S1 < 0 || c_ST < 0) {
            std::fprintf(stderr, "ERROR: %s missing group columns\n", path.c_str());
            std::fclose(fp);
            return 1;
        }
        while (std::fgets(buf, sizeof(buf), fp)) {
            std::vector<std::string> r = split_comma(buf);
            if ((int)r.size() <= c_ST) continue;
            GroupKey k{ r[c_out], r[c_g] };
            double v;
            auto& vals = rows[k];
            vals.S1.push_back(parse_double(r[c_S1], v) ? v : std::nan(""));
            vals.ST.push_back(parse_double(r[c_ST], v) ? v : std::nan(""));
            vals.output_mean.push_back(
                c_om >= 0 && c_om < (int)r.size() && parse_double(r[c_om], v) ? v : std::nan(""));
            vals.output_var.push_back(
                c_ov >= 0 && c_ov < (int)r.size() && parse_double(r[c_ov], v) ? v : std::nan(""));
        }
        std::fclose(fp);
    }

    FILE* fo = std::fopen(out_path, "w");
    if (!fo) { std::fprintf(stderr, "ERROR: cannot write %s\n", out_path); return 1; }
    std::fprintf(fo, "output,group,S1,ST,S1_lo,S1_hi,ST_lo,ST_hi,output_mean,output_variance,n_chunks,S1_se,ST_se\n");
    for (const auto& kv : rows) {
        const auto& k = kv.first;
        const auto& v = kv.second;
        double mS1 = mean_finite(v.S1);
        double mST = mean_finite(v.ST);
        double seS1 = se_finite(v.S1);
        double seST = se_finite(v.ST);
        double mOM = mean_finite(v.output_mean);
        double mOV = mean_finite(v.output_var);
        double S1_lo = mS1 - 1.96*seS1, S1_hi = mS1 + 1.96*seS1;
        double ST_lo = mST - 1.96*seST, ST_hi = mST + 1.96*seST;
        int n_valid = 0;
        for (double x : v.S1) if (std::isfinite(x)) ++n_valid;
        std::fprintf(fo, "%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%.6g,%.6g\n",
            k.output.c_str(), k.group.c_str(),
            mS1, mST, S1_lo, S1_hi, ST_lo, ST_hi,
            mOM, mOV, n_valid, seS1, seST);
    }
    std::fclose(fo);

    std::printf("merged %d chunks into %s (%zu rows)\n",
                n_chunks, out_path, rows.size());
    return 0;
}

void print_usage(const char* prog) {
    std::fprintf(stderr,
        "usage:\n"
        "  %s <out.csv> <chunk0.csv> <chunk1.csv> [chunk2.csv ...]            # first-order\n"
        "  %s --pairs <out_pairs.csv> <c0_pairs.csv> <c1_pairs.csv> [...]     # second-order pairs\n"
        "  %s --groups <out_groups.csv> <c0_groups.csv> <c1_groups.csv> [...] # groups\n",
        prog, prog, prog);
}

}  // anon

int main(int argc, char** argv) {
    if (argc < 4) { print_usage(argv[0]); return 1; }
    enum Mode { FIRST_ORDER, PAIRS, GROUPS } mode = FIRST_ORDER;
    int arg_start = 1;
    if (std::strcmp(argv[1], "--pairs") == 0)  { mode = PAIRS;  arg_start = 2; }
    if (std::strcmp(argv[1], "--groups") == 0) { mode = GROUPS; arg_start = 2; }
    if (argc < arg_start + 3) { print_usage(argv[0]); return 1; }
    const char* out = argv[arg_start];
    std::vector<std::string> chunks;
    for (int i = arg_start + 1; i < argc; ++i) chunks.push_back(argv[i]);
    switch (mode) {
        case FIRST_ORDER: return merge_first_order(out, chunks);
        case PAIRS:       return merge_pairs(out, chunks);
        case GROUPS:      return merge_groups(out, chunks);
    }
    return 1;
}
