import sys
import array
import argparse
import time
import numpy as np

class FMIndex:
    def __init__(self, text, checkpoint_rate=128, sa_sample_rate=32):
        """
        Args:
            text (str): The genome/text to index.
            checkpoint_rate (int): Interval for Occ table checkpoints.
            sa_sample_rate (int): Interval for Suffix Array sampling.
        """
        if not text.endswith('$'):
            text += '$'
        
        self.checkpoint_rate = checkpoint_rate
        self.sa_sample_rate = sa_sample_rate
        
        # 1. Build Suffix Array using SA-IS
        print("  > constructing Suffix Array (SA-IS)...")
        unique_chars = sorted(list(set(text)))
        self.sigma = len(unique_chars)
        self.char_map = {c: i for i, c in enumerate(unique_chars)}
        
        # int_text for SA-IS
        int_text = array.array('I', [self.char_map[c] + 1 for c in text]) 
        sa_list = self._build_sa_is(int_text, self.sigma + 1)
        del int_text


        # 2. Build BWT and SSA
        print("  > building BWT and SSA...")
        n = len(text)
        bwt_list = [''] * n
        self.ssa = {}
        
        for i, sa_val in enumerate(sa_list):
            if sa_val % self.sa_sample_rate == 0:
                self.ssa[i] = sa_val
            
            if sa_val == 0:
                bwt_list[i] = '$'
            else:
                bwt_list[i] = text[sa_val - 1]

        # Store BWT as bytes
        self.bwt_arr = "".join(bwt_list).encode('ascii')
        del bwt_list
        del sa_list

        # 3. Build C Table
        print("  > building C Table...")
        # Count using integers (0-255) directly to avoid chr() overhead
        counts = {}
        for char_code in self.bwt_arr:
            counts[char_code] = counts.get(char_code, 0) + 1
            
        self.c_table = {}
        current = 0
        # Sort by ASCII value
        for char_code in sorted(counts.keys()):
            self.c_table[char_code] = current
            current += counts[char_code]

        # 4. Build Checkpoints (Integer Keys)
        print(f"  > building Occ Checkpoints (Rate: {self.checkpoint_rate})...")
        self.checkpoints = {c: array.array('I') for c in self.c_table.keys()}
        running_counts = {c: 0 for c in self.c_table.keys()}
        
        for i, char_code in enumerate(self.bwt_arr):
            if i % self.checkpoint_rate == 0:
                for c in running_counts:
                    self.checkpoints[c].append(running_counts[c])
            running_counts[char_code] += 1

    def _rank(self, char_code, i):
        if i == 0: return 0
        
        # 1. Get nearest previous checkpoint
        chunk_idx = i // self.checkpoint_rate
        
        # Clamp index to prevent IndexError at the very end of genome
        if chunk_idx >= len(self.checkpoints[char_code]):
            chunk_idx = len(self.checkpoints[char_code]) - 1
            
        checkpoint_val = self.checkpoints[char_code][chunk_idx]
        
        # 2. Scan remaining bytes (Fast C-level count)
        start_idx = chunk_idx * self.checkpoint_rate
        slice_count = self.bwt_arr[start_idx:i].count(char_code)
        
        return checkpoint_val + slice_count

    def count(self, pattern):
        # Convert pattern to bytes ONCE
        pattern_bytes = pattern.encode('ascii')
        
        sp = 0
        ep = len(self.bwt_arr) - 1
        
        # Iterate backwards over integers (bytes)
        for char_code in reversed(pattern_bytes):
            if char_code not in self.c_table:
                return 0, -1
            
            c_val = self.c_table[char_code]
            
            sp = c_val + self._rank(char_code, sp)
            ep = c_val + self._rank(char_code, ep + 1) - 1
            
            if sp > ep:
                return 0, -1
                
        return sp, ep

    def locate(self, pattern):
        sp, ep = self.count(pattern)
        if sp > ep:
            return []
        
        results = []
        for i in range(sp, ep + 1):
            curr_row = i
            steps = 0
            
            # Walk backwards until we hit a sampled SA row
            while curr_row not in self.ssa:
                char_code = self.bwt_arr[curr_row] # int
                
                # LF Step
                curr_row = self.c_table[char_code] + self._rank(char_code, curr_row)
                steps += 1
                
            results.append(self.ssa[curr_row] + steps)
            
        return sorted(results)

    # SA-IS Construction Algorithm
    def _build_sa_is(self, s, K):
        n = len(s)
        t = array.array('b', [0]*n)
        t[n-1] = 1
        for i in range(n-2, -1, -1):
            if s[i] < s[i+1]: t[i] = 1
            elif s[i] > s[i+1]: t[i] = 0
            else: t[i] = t[i+1]
            
        lms_indices = array.array('I')
        is_lms = array.array('b', [0]*n)
        for i in range(1, n):
            if t[i] == 1 and t[i-1] == 0:
                is_lms[i] = 1
                lms_indices.append(i)
                
        def get_buckets(s, K, end=True):
            counts = array.array('I', [0]*K)
            for x in s: counts[x] += 1
            buckets = array.array('I', [0]*K)
            s_val = 0
            for i in range(K):
                s_val += counts[i]
                buckets[i] = s_val if end else s_val - counts[i]
            return buckets

        def induced_sort(s, sa, K, lms_inds=None):
            for i in range(len(sa)): sa[i] = -1
            buckets = get_buckets(s, K, end=True)
            if lms_inds:
                for i in range(len(lms_inds)-1, -1, -1):
                    idx = lms_inds[i]
                    c = s[idx]
                    buckets[c] -= 1
                    sa[buckets[c]] = idx
            
            buckets = get_buckets(s, K, end=False)
            for i in range(n):
                idx = sa[i]
                if idx > 0 and t[idx-1] == 0:
                    c = s[idx-1]
                    sa[buckets[c]] = idx - 1
                    buckets[c] += 1
            
            buckets = get_buckets(s, K, end=True)
            for i in range(n-1, -1, -1):
                idx = sa[i]
                if idx > 0 and t[idx-1] == 1:
                    c = s[idx-1]
                    buckets[c] -= 1
                    sa[buckets[c]] = idx - 1

        sa = array.array('i', [-1]*n)
        induced_sort(s, sa, K, lms_indices)
        
        sorted_lms = array.array('I', [x for x in sa if is_lms[x]])
        lms_names = array.array('i', [-1]*n)
        curr_name = 0
        lms_names[sorted_lms[0]] = 0
        for i in range(1, len(sorted_lms)):
            u = sorted_lms[i-1]
            v = sorted_lms[i]
            diff = False
            for k in range(n):
                if s[u+k] != s[v+k] or t[u+k] != t[v+k]:
                    diff = True; break
                if k > 0 and (is_lms[u+k] or is_lms[v+k]):
                    break
            if diff: curr_name += 1
            lms_names[v] = curr_name
            
        summary_str = array.array('I', [x for x in lms_names if x != -1])
        summary_sa = []
        if curr_name + 1 < len(summary_str):
            summary_sa = self._build_sa_is(summary_str, curr_name + 1)
        else:
            summary_sa = array.array('I', [0] * len(summary_str))
            for i, name in enumerate(summary_str): summary_sa[name] = i
                
        sorted_lms = array.array('I', [lms_indices[x] for x in summary_sa])
        induced_sort(s, sa, K, sorted_lms)
        return sa

# Main & Helper
def read_fasta(filepath):
    full_seq = []
    print(f"Reading {filepath}...")
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('>'): continue
                full_seq.append(line.strip().upper())
        print(f"File read. Joining {len(full_seq)} lines...")
        return "".join(full_seq)
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='Path to FASTA file')
    parser.add_argument('-p', '--pattern', help='Pattern to search', required=True)
    args = parser.parse_args()

    t0 = time.perf_counter()
    genome = read_fasta(args.fasta_file)
    t_read = time.perf_counter() - t0
    
    t0 = time.perf_counter()
    fm = FMIndex(genome)
    t_build = time.perf_counter() - t0

    print(f"\nStats:")
    print(f"  Genome Length: {len(genome)}")
    print(f"  Read Time:     {t_read:.4f}s")
    print(f"  Build Time:    {t_build:.4f}s")
    
    t0 = time.perf_counter()
    sp, ep = fm.count(args.pattern)
    matches = max(0, ep - sp + 1)
    
    locs = []
    if matches > 0 and matches < 50: 
        locs = fm.locate(args.pattern)
    t_query = time.perf_counter() - t0
    
    print(f"\nQuery Results for '{args.pattern}':")
    print(f"  Count: {matches}")
    if locs:
        print(f"  Locations: {locs}")
    elif matches >= 50:
        print(f"  Locations: (Too many to list, found {matches})")
    print(f"  Time:  {t_query:.2e}s")