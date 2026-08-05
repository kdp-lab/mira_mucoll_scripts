import re, sys, argparse
from datetime import datetime
from collections import defaultdict

TS = "%Y-%m-%d %H:%M:%S"
PAT = re.compile(r'^ *040 \((\d+)\.(\d+)\.(\d+)\) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}).*Finished transferring (input|output) files')

def hms(sec):
    s=int(round(sec)); h=s//3600; m=(s%3600)//60; return f"{h}h {m:02d}m {s%60:02d}s" if h else (f"{m}m {s%60:02d}s" if m else f"{s}s")

def stream(files):
    if not files:
        for line in sys.stdin: yield line
    else:
        for f in files:
            with open(f, 'r', encoding='utf-8', errors='replace') as fh:
                for line in fh: yield line

def main():
    ap = argparse.ArgumentParser(description="Bin HTCondor IO transfer durations by proc ranges.")
    ap.add_argument("logs", nargs="*")
    ap.add_argument("--bin", type=int, default=100, help="proc bin size (default 100)")
    ap.add_argument("--csv", action="store_true", help="CSV output (one row per bin)")
    args = ap.parse_args()

    pending = {}  
    bins = defaultdict(list) 

    for line in stream(args.logs):
        m = PAT.match(line)
        if not m: 
            continue
        cid, proc, sub, ts_str, io = m.groups()
        attempt = f"{cid}.{proc}.{sub}"
        t = datetime.strptime(ts_str, TS)

        if io == "input":
            pending[attempt] = t
        else:
            start = pending.pop(attempt, None)
            if not start or t < start:
                continue
            sec = (t - start).total_seconds()

            p = int(proc)
            b0 = (p // args.bin) * args.bin
            b1 = b0 + args.bin - 1
            key = (cid, b0, b1)  
            bins[key].append(sec)

    if not bins:
        print("No input/output transfer pairs found.")
        return

    if args.csv:
        print("cluster,proc_bin_start,proc_bin_end,count,avg_s,min_s,max_s")
        for (cid, b0, b1) in sorted(bins, key=lambda k: (int(k[0]), k[1])):
            vals = bins[(cid,b0,b1)]
            print(f"{cid},{b0},{b1},{len(vals)},{sum(vals)/len(vals):.3f},{min(vals):.3f},{max(vals):.3f}")
    else:
        for (cid, b0, b1) in sorted(bins, key=lambda k: (int(k[0]), k[1])):
            vals = bins[(cid,b0,b1)]
            avg = sum(vals)/len(vals)
            print(f"cluster {cid}  proc {b0:>4d}-{b1:<4d}  "
                  f"n={len(vals):3d}  avg={hms(avg)}  "
                  f"min={hms(min(vals))}  max={hms(max(vals))}")

if __name__ == "__main__":
    main()
