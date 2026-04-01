from collections import namedtuple, defaultdict, Counter
import gzip
import pandas as pd
from tqdm import tqdm
import numpy as np 
from matplotlib import pyplot as plt
import zipfile
import io
from fast_edit_distance import edit_distance, sub_edit_distance
from io import StringIO
import multiprocessing as mp
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import networkx as nx

# a light class for a read in fastq file
read_tuple = namedtuple('read_tuple', ['id', 'seq', 'q_letter'])
def fastq_parser(file_handle):
    while True:
        id = next(file_handle, None)
        if id is None:
            break
        seq = next(file_handle)
        next(file_handle) # skip  '+'
        q_letter = next(file_handle)
        yield read_tuple(id[1:].split()[0], seq.strip(), q_letter.strip()) #每次yield一条read的信息

# split any iterator in to batches  
def batch_iterator(iterator, batch_size):
    """generateor of batches of items in a iterator with batch_size.
    """
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(entry)
        
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):  #保证批次处理的时候，最后一批不满足batch_size的那些数据，也可以被yield
        yield batch

def read_batch_generator(fastq_fns, batch_size):   #输出batch size read info
    """Generator of barches of reads from list of fastq files

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq = fastq_parser(handle)
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            with open(fn, "r") as handle:
                fastq = fastq_parser(handle)
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
def reverse_complement(seq):
    '''
    Args: <str>
        queried seq
    Returns: <str>
        reverse_complement seq
    '''
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                    'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    letters = \
        [comp[base] if base in comp.keys() else base for base in seq]
    return ''.join(letters)[::-1]

def polyA_trimming_idx(seq, seed="AAAA", window=10, min_A=7, min_tail_len=8):
    """
    从 read 末端往前检测 polyA，返回 polyA 起始的绝对坐标（0-based）。
    若未检测到则返回 None。
    """
    s = seq.upper()
    anchor = s.rfind(seed)  # 最右侧 seed 的起点（绝对坐标）
    if anchor == -1:
        return None

    polyA_start = anchor  # 先把 polyA 起点放在 seed 起点
    i = anchor - 1        # 从 seed 之前的碱基开始，向左延伸

    while i >= 0:
        if s[i] == 'A':
            polyA_start = i
            i -= 1
            continue
        left = max(0, i - window + 1)
        if s[left:i+1].count('A') >= min_A:
            polyA_start = i
            i -= 1
            continue
        break

    # 确保 polyA 足够长
    if len(s) - polyA_start < min_tail_len:
        return None
    return polyA_start

def polyA_trimming_idx_neg(seq, **kwargs):
    idx_abs = polyA_trimming_idx(seq, **kwargs)  # 用上面的绝对坐标函数
    if idx_abs is None:
        return None
    return idx_abs - len(seq)  # 负数：从末尾往前的偏移

def find_pos(s, sub):
    return s.find(sub)         # 找不到就是 -1


def rfind_with_negative(s, sub):
    pos = s.rfind(sub)
    if pos == -1:
        return -1  # not found 
    return pos - len(s)

def default_count_threshold_calculation(count_array, exp_cells):
    top_count = np.sort(count_array)[::-1][:exp_cells]
    return np.quantile(top_count, 0.95)/20

def knee_plot(counts, threshold=None, out_fn = 'knee_plot.png'):
    """
    Plot knee plot using the high-confidence putative BC counts

    Args:
        counts (list): high-confidence putative BC counts
        threshold (int, optional): a line to show the count threshold. Defaults to None.
    """
    counts = sorted(counts)[::-1]
    plt.figure(figsize=(8, 8))
    plt.title(f'Barcode rank plot (from high-quality putative BC)')
    plt.loglog(counts,marker = 'o', linestyle="", alpha = 1, markersize=6)
    plt.xlabel('Barcodes')
    plt.ylabel('Read counts')
    plt.axhline(y=threshold, color='r', linestyle='--', label = 'cell calling threshold')
    plt.legend()
    plt.savefig(out_fn)

def get_bc_whitelist(raw_bc_count, full_bc_whitelist=None, exp_cells=None, out_plot_fn = None,empty_max_count = np.inf, DEFAULT_EMPTY_DROP_MIN_ED=None, DEFAULT_EMPTY_DROP_NUM=None):
    percentile_count_thres = default_count_threshold_calculation
    whole_whitelist = []
    if full_bc_whitelist.endswith('.zip'):
        with zipfile.ZipFile(full_bc_whitelist) as zf:
            # check if there is only 1 file
            assert len(zf.namelist()) == 1

            with io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding="utf-8") as f:
                for line in f:
                    whole_whitelist.append(line.strip())
    elif full_bc_whitelist.endswith(".gz"):
        with gzip.open(full_bc_whitelist, "rt", encoding="utf-8") as f:
            for line in f:
                whole_whitelist.append(line.strip())
    else:
        with open(full_bc_whitelist, 'r') as f:
            for line in f:
                whole_whitelist.append(line.strip())
    
    whole_whitelist = set(whole_whitelist)
    raw_bc_count = {k:v for k,v in  raw_bc_count.items() if k in whole_whitelist}
    #print(len(raw_bc_count))`
    t = percentile_count_thres(list(raw_bc_count.values()), exp_cells) #t是
    knee_plot(list(raw_bc_count.values()), t, out_plot_fn)
    cells_bc = {k:v for k,v in raw_bc_count.items() if v > t}
    
    ept_bc = []
    ept_bc_max_count = min(cells_bc.values()) #空bc最大的read支持数
    ept_bc_max_count = min(ept_bc_max_count, empty_max_count)
    #print(ept_bc_max_count)

    ept_bc_candidate = [k for k,v in raw_bc_count.items() if v < ept_bc_max_count]
    #print(len(ept_bc_candidate))
    for k in ept_bc_candidate:
        if min([edit_distance(k, x, max_ed = DEFAULT_EMPTY_DROP_MIN_ED) for x in cells_bc.keys()]) >= DEFAULT_EMPTY_DROP_MIN_ED:
            ept_bc.append(k)
            #print(len(ept_bc))
        # we don't need too much BC in this list
        if len(ept_bc) >  DEFAULT_EMPTY_DROP_NUM:
            break
    return cells_bc, ept_bc

class read_fastq:
    """This class is for mimic the Bio.SeqIO fastq record. The SeqIO is not directly used because it's slow.
    """
    def __init__(self, title, sequence, qscore, quality_map = False):
        self.id = title.split()[0].strip("@")
        self.seq = sequence
        self.qscore = qscore

def _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batch_size):
    """Generator of barches of reads from list of fastq files with the idx of the first read
    in each batch

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 1000.
    """
    read_idx = 0
    putative_bc_f = open(putative_bc_csv, 'r')
    putative_bc_header = next(putative_bc_f)

    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in fastq_parser(handle))

                batch_iter = batch_iterator(fastq, batch_size=batch_size)
                
                for batch in batch_iter:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO(
                            putative_bc_header + \
                            ''.join([next(putative_bc_f) for x in range(batch_len)])
                        ))
                    yield batch, read_idx, batch_bc_df
                    read_idx += batch_len
        else:
            with open(fn) as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in fastq_parser(handle))
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO( putative_bc_header + \
                        ''.join([next(putative_bc_f) for x in range(batch_len)]))
                    )

                    yield batch, read_idx, batch_bc_df #[read_id_seq_qv,....,] 0,1000,2000... df
                    read_idx += batch_len
    putative_bc_f.close()

def _match_bc_row_dual(row, wl3, wl5, max_ed, minQ=0):
    strand = '+'

    # 取出两端字段（你 df 里列名按你说的）
    bc3 = getattr(row, "putative_bc", "") or ""
    umi3 = getattr(row, "putative_umi", "") or ""
    bc5 = getattr(row, "putative_bc_5p", "") or ""
    umi5 = getattr(row, "putative_umi_5p", "") or ""

    # 可选：质量过滤（你如果有两端 qscore，就分别判断；没有就忽略）
    # q3 = getattr(row, "putative_bc_qscore", None)
    # q5 = getattr(row, "putative_bc_qscore_5p", None)
    # if minQ and ((q3 is not None and q3 < minQ) and (q5 is not None and q5 < minQ)):
    #     return ['', '', strand, '', '']

    # 只把“在白名单”当 ok（更符合语义）
    ok3 = (bc3 in wl3)
    ok5 = (bc5 in wl5)

    # 3 ok, 5 不 ok：纠错5
    if ok3 and (not ok5):
        bc5_corr = correct_bc_one_side(bc5, wl5, max_ed)
        umi5_out = umi5 if bc5_corr != '' else ''   # 纠错失败就清空 umi5
        return [bc3, umi3, strand, bc5_corr, umi5_out]

    # 5 ok, 3 不 ok：纠错3
    if ok5 and (not ok3):
        bc3_corr = correct_bc_one_side(bc3, wl3, max_ed)
        umi3_out = umi3 if bc3_corr != '' else ''   # 纠错失败就清空 umi3
        return [bc3_corr, umi3_out, strand, bc5, umi5]

    # 两端都不 ok：两端都纠错
    bc3_corr = correct_bc_one_side(bc3, wl3, max_ed)
    bc5_corr = correct_bc_one_side(bc5, wl5, max_ed)
    umi3_out = umi3 if bc3_corr != '' else ''
    umi5_out = umi5 if bc5_corr != '' else ''
    return [bc3_corr, umi3_out, strand, bc5_corr, umi5_out]

def assign_read_batches(r_batch,
                        whitelist_3p, whitelist_5p,
                        max_ed, gz, minQ=0,
                        emit_unmatched_fastq=True):

    """
    whitelist_3p:指的是whitelist_3p.csv
    whitelist_5p:指的是whitelist_5p.csv
    """
    read_batch, start_df_idx, df = r_batch
    df = df.fillna('')

    wl3 = set(whitelist_3p)
    wl5 = set(whitelist_5p)

    out_buffer = ''
    unmatched_fastq_buffer = ''

    # 1) 对每行同时纠错 3' 和 5'，生成 5 列
    new_cols = []
    for row in df.itertuples():
        new_cols.append(_match_bc_row_dual(row, wl3, wl5, max_ed, minQ))

    df[['BC3_corrected', 'putative_umi', 'strand',
        'BC5_corrected', 'putative_umi_5p']] = new_cols #在之前的putative上新增了这五列

    # 2) 统计“成功”的 read（你可以定义为：至少一端成功）
    ok3 = (df['BC3_corrected'] != '') & (df['putative_umi'] != '') #ok3和ok5是布尔值
    ok5 = (df['BC5_corrected'] != '') & (df['putative_umi_5p'] != '')
    demul_read_count = int((ok3 | ok5).sum()) #只有一端拆到就算成功的数量，比例应该不会大于1

    # 3) 写 fastq：至少一端 (BC+UMI) 成功才输出
    for r, bc in zip(read_batch, df.itertuples()):
        try:
            assert bc.read_id == r.id
        except AssertionError:
            err_msg("Different order in putative bc file and input fastq!", printit=True)
            sys.exit()

        side3_ok = (bc.BC3_corrected != '' and bc.putative_umi != '') #putative.csv中3‘端是否成功
        side5_ok = (bc.BC5_corrected != '' and bc.putative_umi_5p != '') #putative.csv中5‘端是否成功

        if not (side3_ok or side5_ok):
            # unmatched（沿用你原来的 A比例过滤逻辑）
            if emit_unmatched_fastq:
                putative_bc = getattr(bc, "putative_bc", "")
                if putative_bc:
                    a_ratio = putative_bc.count("A") / len(putative_bc)
                    if a_ratio <= 0.5:
                        cb3 = bc.BC3_corrected if bc.BC3_corrected else "NA"
                        ub3 = bc.putative_umi if bc.putative_umi else "NA"
                        cb5 = bc.BC5_corrected if bc.BC5_corrected else "NA"
                        ub5 = bc.putative_umi_5p if bc.putative_umi_5p else "NA"
                        #header = (f"@{cb3}_{ub3}|{cb5}_{ub5}#{bc.read_id}_{getattr(bc,'strand','+')}"
                                  #f"\tCB3:Z:{cb3}\tUB3:Z:{ub3}\tCB5:Z:{cb5}\tUB5:Z:{ub5}")
                        #header = (
                        #    f"@{bc.read_id}"
                        #    f"\tCB:Z:NA\tUB:Z:NA"
                        #    f"\tCB3:Z:{cb3}\tUB3:Z:{ub3}\tCB5:Z:{cb5}\tUB5:Z:{ub5}"
                        #    f"\tSTR:Z:{getattr(bc,'strand','+')}"
                        #)
                        header = f"@{bc.read_id}\n"
                        unmatched_fastq_buffer += header
                        unmatched_fastq_buffer += r.seq + '\n+\n' + r.qscore + '\n'
            continue

        # 4) 裁剪策略：优先用 3' 的 polyA/umi_fixed；没有就退到 5' umi_fixed（如果你有该列）
        
        has3 = (getattr(bc, "BC3_corrected", "") != "")
        has5 = (getattr(bc, "BC5_corrected", "") != "")
        
        # ======= 更稳健：只有 barcode+UMI 都有，才允许裁剪该端 =======
        trim3_ok = has3 and (getattr(bc, "putative_umi", "") != "")
        trim5_ok = has5 and (getattr(bc, "putative_umi_5p", "") != "")
        
        seq = r.seq
        qscore = r.qscore
        L = len(seq)

        umi5_start = _to_int(getattr(bc, "umi_fixed_locs_5p", None))   # 5' UMI 起点
        umi5 = getattr(bc, "putative_umi_5p", "") or ""
        start_cut = 0
        if trim5_ok and umi5_start is not None:
            start_cut = umi5_start + len(umi5)
        
        start_cut = max(0, min(L, start_cut))

        polyA = _to_int(getattr(bc, "polyA_starts", None))
        umi3_loc = _to_int(getattr(bc, "umi_fixed_locs", None))
        end_cut = L
        if trim3_ok:
            # 1) 优先 polyA（如果是负数，Python 里 end_cut = L + polyA）
            if polyA is not None:
                end_cut = (L + polyA) if polyA < 0 else polyA
            # 2) 没 polyA 用 umi_fixed_locs（沿用你原经验：cut = umi_fixed_locs - 10）
            elif umi3_loc is not None:
                cut = umi3_loc - 10
                end_cut = (L + cut) if cut < 0 else cut
            else:
                end_cut = L
        end_cut = max(0, min(L, end_cut))

        if start_cut >= end_cut:
            continue  # 剪坏/剪空了就放弃（也可以写到 unmatched 里）
        
        seq = seq[start_cut:end_cut]
        qscore = qscore[start_cut:end_cut]
        
        # 可选：太短的不输出
        if len(seq) < 30:
            continue

        # 5) header：把两端都写进去；同时给一个“主 CB/UB”（优先 3'，否则用 5'）
        """
        if side3_ok:
            CB = bc.BC3_corrected
            UB = bc.putative_umi
        else:
            CB = bc.BC5_corrected
            UB = bc.putative_umi_5p
        """

        cb3 = bc.BC3_corrected if bc.BC3_corrected else "NA"
        ub3 = bc.putative_umi if bc.putative_umi else "NA"
        cb5 = bc.BC5_corrected if bc.BC5_corrected else "NA"
        ub5 = bc.putative_umi_5p if bc.putative_umi_5p else "NA"

        if side3_ok:
            cb_main = cb3
            umi_main = ub3
        else:
            cb_main = cb5
            umi_main = ub5

        umi_main = umi_main.replace(" ", "").upper()
        if any(ch not in "ACGTN" for ch in umi_main):
            continue

        #read_name = f"{bc.read_id}_{umi_main}"
        read_name = bc.read_id

        # 你想把它们“全部输出在 readid 上”
        #read_name = f"{cb3}_{ub3}|{cb5}_{ub5}#{bc.read_id}_{bc.strand}"
        
        #out_buffer += (
        #    f"@{read_name}"
        #    f"\tCB:Z:{cb_main}\tUB:Z:{umi_main}"
        #    f"\tCB3:Z:{cb3}\tUB3:Z:{ub3}\tCB5:Z:{cb5}\tUB5:Z:{ub5}"
        #    f"\tSTR:Z:{bc.strand}\n"
        #)
        out_buffer += f"@{read_name}\n"
        out_buffer += seq + "\n+\n" + qscore + "\n"

    # 6) gzip / plain
    b_out_buffer = gzip.compress(out_buffer.encode('utf-8')) if gz else out_buffer.encode('utf-8')
    if emit_unmatched_fastq:
        b_unmatched_fastq = gzip.compress(unmatched_fastq_buffer.encode('utf-8')) if gz else unmatched_fastq_buffer.encode('utf-8')
    else:
        b_unmatched_fastq = None

    return df, b_out_buffer, demul_read_count, len(read_batch), b_unmatched_fastq

def _to_int(x):
    if x in ('', None):
        return None
    try:
        return int(float(x))  # 兼容 26.0
    except Exception:
        return None

def err_msg(msg, printit = False):
    CRED = '\033[91m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def warning_msg(msg, printit = False):
    CRED = '\033[93m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def green_msg(msg, printit = False):
    CRED = '\033[92m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND


def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,
                           pbar=True, pbar_unit='Read',pbar_func=len, 
                           schduler = 'process', *arg, **kwargs):
    """multiple processing or threading, 

    Args:
        func: function to be run parallely
        iterator: input to the function in each process/thread
        n_process (int, optional): number of cores or threads. Defaults to mp.cpu_count()-1.
        pbar (bool, optional): Whether or not to output a progres bar. Defaults to True.
        pbar_unit (str, optional): Unit shown on the progress bar. Defaults to 'Read'.
        pbar_func (function, optional): Function to calculate the total length of the progress bar. Defaults to len.
        schduler (str, optional): 'process' or 'thread'. Defaults to 'process'.

    Yields:
        return type of the func: the yield the result in the order of submit
    """
    class fake_future:
        # a fake future class to be used in single processing
        def __init__(self, rst):
            self.rst = rst
        def result(self):
            return self.rst

    if schduler == 'process':
        # make sure the number of process is not larger than the number of cores
        n_process = min(n_process-1, mp.cpu_count()-1)
        if n_process > 1:
            executor = concurrent.futures.ProcessPoolExecutor(n_process)
    elif schduler == 'thread':
        if n_process > 1:
            executor = concurrent.futures.ThreadPoolExecutor(n_process)
    else:
        green_msg('Error in multiprocessing_submit: schduler should be either process or thread', printit=True)
        sys.exit(1)

    if pbar:
        _pbar = tqdm(unit=pbar_unit, desc='Processed')
        
    # run in single process/thread if n_process < 1
    if n_process <= 1:
        for it in iterator:
            yield fake_future(func(it, *arg, **kwargs))
            if pbar:
                _pbar.update(pbar_func(it))
        return

    # A dictionary which will contain the future object
    max_queue = n_process
    futures = {}
    n_job_in_queue = 0
    
    # make sure the result is yield in the order of submit.
    job_idx = 0
    job_completed = {}

    # submit the first batch of jobs
    while n_job_in_queue < max_queue:
        i = next(iterator, None)
        if i is None:
            break
        futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
        job_idx += 1
        n_job_in_queue += 1
        job_to_yield = 0
    # yield the result in the order of submit and submit new jobs
    while True:
        # will wait until as least one job finished
        # batch size as value, release the cpu as soon as one job finished
        job = next(as_completed(futures), None)

        # yield the completed job in the order of submit  
        if job is not None:
            job_completed[futures[job][1]] = job, futures[job][0]
            del futures[job]

        # 
        if job is None and i is None and len(job_completed)==0:
            break

        # check order
        while job_to_yield in job_completed.keys():
            # update pregress bar based on batch size
            if pbar:
                _pbar.update(job_completed[job_to_yield][1])
            yield job_completed[job_to_yield][0]
            del job_completed[job_to_yield]
            
            # submit new job
            i = next(iterator, None)
            if i is not None:
                futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
                job_idx += 1
                
            job_to_yield += 1

def correct_bc_one_side(bc, whitelist, max_ed):
    """return corrected_bc or '' (failed); if bc is '' -> ''"""
    if not bc:  # 空字符串 / None
        return ''
    if bc in whitelist:
        return bc

    best_ed = max_ed
    bc_hit = ''

    for i in whitelist:
        ed, _ = sub_edit_distance(i, bc, best_ed)
        if ed < best_ed:
            best_ed = ed
            bc_hit = i
        elif ed == best_ed:
            if not bc_hit:
                bc_hit = i
            else:
                bc_hit = 'ambiguous'
                best_ed -= 1
                if best_ed < 0:
                    return ''
    if bc_hit in ('', 'ambiguous'):
        return ''
    return bc_hit


def assign_read(fastq_fns=None, fastq_out=None, putative_bc_csv=None, 
                    whitelsit_3p=None, whitelsit_5p=None, max_ed=None,n_process=None, batchsize=None, minQ=0):
    
    gz = fastq_out.endswith('.gz') #判断输出文件是否为gz文件
    out_dir = os.path.dirname(fastq_out)       # 输出目录
    unmatched_out = os.path.join(out_dir, "unmatched_reads.fastq.gz")
        
    r_batches = \
        _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batchsize)
    
    whitelist_3p_list = [] 
    with open(whitelsit_3p, 'r') as f:
        for line in f:
            whitelist_3p_list.append(line.split('-')[0].strip())

    whitelist_5p_list = [] 
    with open(whitelsit_5p, 'r') as f:
        for line in f:
            whitelist_5p_list.append(line.split('-')[0].strip())

    if n_process == 1:
        demul_count_tot = 0
        count_tot = 0
        with open(fastq_out, 'wb') as output_handle, open(unmatched_out, 'wb') as unmatched_handle:
            pbar = tqdm(unit="Reads", desc='Processed')
            for r_batch in r_batches:
                _, b_fast_str, demul_count, read_count, b_unmatched_fastq = assign_read_batches(r_batch, whitelist_3p_list, whitelist_5p_list, max_ed,  gz, minQ=0, emit_unmatched_fastq=True)
                demul_count_tot += demul_count
                count_tot += read_count
                output_handle.write(b_fast_str)
                if b_unmatched_fastq:
                    unmatched_handle.write(b_unmatched_fastq)
                pbar.update(read_count) #
        green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!")
        green_msg(f"Unmatched reads saved in: {unmatched_out}!")
        
    else:
        rst_futures = multiprocessing_submit(assign_read_batches, 
                            r_batches, 
                            n_process=n_process,
                            schduler = "process",
                            pbar_func=lambda x: len(x[0]),
                            whitelist_3p = whitelist_3p_list,
                            whitelist_5p = whitelist_5p_list,
                            max_ed = max_ed,
                            gz = gz)

        demul_count_tot = 0
        count_tot = 0
        df_list = []
        with open(fastq_out, 'wb') as output_handle, open(unmatched_out, 'wb') as unmatched_handle:
            for f in rst_futures:
                df, b_fast_str, demul_count, read_count, b_unmatched_fastq = f.result()
                demul_count_tot += demul_count
                count_tot += read_count
                output_handle.write(b_fast_str)

                if b_unmatched_fastq:
                    unmatched_handle.write(b_unmatched_fastq)
                # 保存df
                df_list.append(df)
        big_df = pd.concat(df_list, ignore_index=True)
    
        green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!")
        green_msg(f"Unmatched reads saved in: {unmatched_out}!")
    
    return demul_count_tot, count_tot,big_df


def _safe_min_q(qstr):
    if qstr is None or len(qstr) == 0:
        return None
    return min([ord(x) for x in qstr]) - 33

def get_3p_features_new(
    read_info,
    read_ids,
    putative_bcs,
    bc_fixed_locs,
    putative_bc_min_qs,
    umis,
    umi_fixed_locs,
    post_umi_flankings,
    polyA_starts,
    read_types,
    FIX6_A="GCTACC",
    FIX5_B="AGATC",
    FIX6_UMI="TAGGCT",
):
    """
    3p structure:
    UMI(10) - FIX6_UMI(6) - BC1(10) - FIX5_B(5) - BC2(10) - FIX6_A(6) - BC3(10)

    putative_bcs: 41bp = BC1 + FIX5_B + BC2 + FIX6_A + BC3
    umis: 16bp = UMI(10) + FIX6_UMI(6)
    """
    part_id = read_info.id
    part_seq = read_info.seq[-60:]          # 新结构更长，窗口放大
    part_qv = read_info.q_letter[-60:]

    read_ids.append(part_id)

    putative_bc_min_q = None
    umi = None
    umi_fixed_loc = None
    post_umi_flanking = None
    polyA_start = None
    read_type = None

    # 理论位置：末端 BC3 为 10bp，因此 FIX6_A 应在 -16
    BC_fixed_loc = rfind_with_negative(part_seq, FIX6_A)
    bc_fixed_locs.append(BC_fixed_loc)

    def _append_polyA_and_type(polyA_anchor_loc, hit_type, miss_type):
        nonlocal polyA_start, read_type
        seq_polyA = read_info.seq[polyA_anchor_loc - 100: polyA_anchor_loc]
        last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
        if last_polyA_idx:
            polyA_start = last_polyA_idx + polyA_anchor_loc
            polyA_starts.append(polyA_start)
            read_type = hit_type
            read_types.append(read_type)
        else:
            polyA_starts.append(polyA_start)
            read_type = miss_type
            read_types.append(read_type)

    if BC_fixed_loc == -16:
        # 完整命中外侧 6bp 固定序列
        fix5_expected = read_info.seq[-31:-26]

        if fix5_expected == FIX5_B:
            putative_bc = read_info.seq[-41:]
            putative_bcs.append(putative_bc)
            putative_bc_min_q = _safe_min_q(read_info.q_letter[-41:])
            putative_bc_min_qs.append(putative_bc_min_q)

            # 在 BC1 前的局部窗口找 UMI 前 6bp 固定序列
            # BC1 起点 = -41
            find_umi_seq = read_info.seq[-57:-41]
            umi_fixed_loc_re = rfind_with_negative(find_umi_seq, FIX6_UMI)

            if umi_fixed_loc_re != -1:
                umi_fixed_loc = umi_fixed_loc_re - 41
                umi_fixed_locs.append(umi_fixed_loc)

                umi = read_info.seq[umi_fixed_loc - 10: umi_fixed_loc + 6]
                umis.append(umi)

                post_umi_flanking = read_info.seq[umi_fixed_loc - 15: umi_fixed_loc - 10]
                post_umi_flankings.append(post_umi_flanking)

                _append_polyA_and_type(umi_fixed_loc - 10, 1, 2)
            else:
                # 找不到 UMI 固定序列，粗略截取 16bp
                umi = read_info.seq[-57:-41]
                umis.append(umi)
                umi_fixed_locs.append(umi_fixed_loc)

                post_umi_flanking = read_info.seq[-62:-57]
                post_umi_flankings.append(post_umi_flanking)

                _append_polyA_and_type(-41, 3, 4)
        else:
            # 找到 FIX6_A 但 FIX5_B 不匹配，仍输出末端 41bp 作为候选
            putative_bc = read_info.seq[-41:]
            putative_bcs.append(putative_bc)
            putative_bc_min_q = _safe_min_q(read_info.q_letter[-41:])
            putative_bc_min_qs.append(putative_bc_min_q)

            find_umi_seq = read_info.seq[-57:-41]
            umi_fixed_loc_re = rfind_with_negative(find_umi_seq, FIX6_UMI)

            if umi_fixed_loc_re != -1:
                umi_fixed_loc = umi_fixed_loc_re - 41
                umi_fixed_locs.append(umi_fixed_loc)

                umi = read_info.seq[umi_fixed_loc - 10: umi_fixed_loc + 6]
                umis.append(umi)

                post_umi_flanking = read_info.seq[umi_fixed_loc - 15: umi_fixed_loc - 10]
                post_umi_flankings.append(post_umi_flanking)

                _append_polyA_and_type(umi_fixed_loc - 10, 5, 6)
            else:
                umi = read_info.seq[-57:-41]
                umis.append(umi)
                umi_fixed_locs.append(umi_fixed_loc)

                post_umi_flanking = read_info.seq[-62:-57]
                post_umi_flankings.append(post_umi_flanking)

                _append_polyA_and_type(-41, 7, 8)

    elif BC_fixed_loc < -16:
        # FIX6_A 更靠左，说明整段 barcode 也更靠左，但大概率仍完整
        bc_start = BC_fixed_loc - 25         # BC1 start
        bc_end = BC_fixed_loc + 16           # BC3 end
        fix5_start = BC_fixed_loc - 15       # FIX5_B start

        putative_bc = read_info.seq[bc_start:bc_end]
        putative_bcs.append(putative_bc)
        putative_bc_min_q = _safe_min_q(read_info.q_letter[bc_start:bc_end])
        putative_bc_min_qs.append(putative_bc_min_q)

        fix5_seq = read_info.seq[fix5_start:fix5_start + 5]

        find_umi_seq = read_info.seq[bc_start - 16:bc_start]
        umi_fixed_loc_re = rfind_with_negative(find_umi_seq, FIX6_UMI)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + bc_start - 16
            umi_fixed_locs.append(umi_fixed_loc)

            umi = read_info.seq[umi_fixed_loc - 10: umi_fixed_loc + 6]
            umis.append(umi)

            post_umi_flanking = read_info.seq[umi_fixed_loc - 15: umi_fixed_loc - 10]
            post_umi_flankings.append(post_umi_flanking)

            if fix5_seq == FIX5_B:
                _append_polyA_and_type(umi_fixed_loc - 10, 9, 10)
            else:
                _append_polyA_and_type(umi_fixed_loc - 10, 11, 12)
        else:
            umi = read_info.seq[bc_start - 16:bc_start]
            umis.append(umi)
            umi_fixed_locs.append(umi_fixed_loc)

            post_umi_flanking = read_info.seq[bc_start - 21:bc_start - 16]
            post_umi_flankings.append(post_umi_flanking)

            if fix5_seq == FIX5_B:
                _append_polyA_and_type(bc_start, 13, 14)
            else:
                _append_polyA_and_type(bc_start, 15, 16)

    elif BC_fixed_loc == -1:
        # 完全找不到外侧 FIX6_A，末端粗略输出 41bp，同时单独找 umi fixed
        putative_bcs.append(read_info.seq[-41:])
        putative_bc_min_qs.append(_safe_min_q(read_info.q_letter[-41:]))

        find_umi_seq = read_info.seq[:-41]
        umi_fixed_loc_re = rfind_with_negative(find_umi_seq, FIX6_UMI)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re
            umi_fixed_locs.append(umi_fixed_loc)

            umi = read_info.seq[umi_fixed_loc - 10: umi_fixed_loc + 6]
            umis.append(umi)

            post_umi_flanking = read_info.seq[umi_fixed_loc - 15: umi_fixed_loc - 10]
            post_umi_flankings.append(post_umi_flanking)

            _append_polyA_and_type(umi_fixed_loc - 10, 17, 18)
        else:
            umi_fixed_locs.append(umi_fixed_loc)
            umis.append(umi)
            post_umi_flankings.append(post_umi_flanking)
            polyA_starts.append(polyA_start)

            read_type = 19
            read_types.append(read_type)

    else:
        # BC_fixed_loc > -16，说明 BC3 不完整，尝试依靠 FIX5_B / UMI 固定序列保守输出
        bc_start = BC_fixed_loc - 25
        fix5_start = BC_fixed_loc - 15
        fix5_seq = read_info.seq[fix5_start:fix5_start + 5]

        find_umi_seq = read_info.seq[bc_start - 16:bc_start]
        umi_fixed_loc_re = rfind_with_negative(find_umi_seq, FIX6_UMI)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + bc_start - 16
            umi_fixed_locs.append(umi_fixed_loc)

            umi = read_info.seq[umi_fixed_loc - 10: umi_fixed_loc + 6]
            umis.append(umi)

            post_umi_flanking = read_info.seq[umi_fixed_loc - 15: umi_fixed_loc - 10]
            post_umi_flankings.append(post_umi_flanking)

            # barcode 不完整时，从 BC1 推到 read 末端
            putative_bc_start_loc = umi_fixed_loc + 6
            putative_bc = read_info.seq[putative_bc_start_loc:]
            putative_bcs.append(putative_bc)
            putative_bc_min_q = _safe_min_q(read_info.q_letter[putative_bc_start_loc:])
            putative_bc_min_qs.append(putative_bc_min_q)

            if fix5_seq == FIX5_B:
                _append_polyA_and_type(umi_fixed_loc - 10, 20, 21)
            else:
                _append_polyA_and_type(umi_fixed_loc - 10, 22, 23)
        else:
            putative_bcs.append(read_info.seq[-41:])
            putative_bc_min_qs.append(_safe_min_q(read_info.q_letter[-41:]))

            umi_fixed_locs.append(umi_fixed_loc)
            umis.append(umi)
            post_umi_flankings.append(post_umi_flanking)
            polyA_starts.append(polyA_start)

            read_type = 24
            read_types.append(read_type)

    return (
        read_ids,
        putative_bcs,
        bc_fixed_locs,
        putative_bc_min_qs,
        umis,
        umi_fixed_locs,
        post_umi_flankings,
        polyA_starts,
        read_types,
    )


def get_5p_features_new(
    read_info,
    read_ids_5p,
    putative_bcs_5p,
    bc_fixed_locs_5p,
    putative_bc_min_qs_5p,
    umis_5p,
    umi_fixed_locs_5p,
    FIX6_A_5p="CCTTCC",
    FIX5_B_5p="TGCTG",
    FIX6_UMI_5p="CCACTG",
):
    """
    5p structure:
    BC1(10) - FIX6_A(6) - BC2(10) - FIX5_B(5) - BC3(10) - FIX6_UMI(6) - UMI(10)

    putative_bcs_5p: 41bp = BC1 + FIX6_A + BC2 + FIX5_B + BC3
    umis_5p: 16bp = FIX6_UMI(6) + UMI(10)
    """
    part_id = read_info.id
    part_seq = read_info.seq[:60]
    part_qv = read_info.q_letter[:60]

    read_ids_5p.append(part_id)

    putative_bc_min_q = None
    umi = None
    umi_fixed_loc = None

    # 理论位置：FIX6_A_5p 应在第 10 位开始
    BC_fixed_loc = find_pos(part_seq, FIX6_A_5p)
    bc_fixed_locs_5p.append(BC_fixed_loc)

    if BC_fixed_loc == 10:
        fix5_expected = read_info.seq[26:31]

        if fix5_expected == FIX5_B_5p:
            putative_bc = read_info.seq[:41]
            putative_bcs_5p.append(putative_bc)
            putative_bc_min_q = _safe_min_q(read_info.q_letter[:41])
            putative_bc_min_qs_5p.append(putative_bc_min_q)

            find_umi_seq = read_info.seq[41:57]
            umi_fixed_loc_re = find_pos(find_umi_seq, FIX6_UMI_5p)

            if umi_fixed_loc_re != -1:
                umi_fixed_loc = umi_fixed_loc_re + 41
                umi_fixed_locs_5p.append(umi_fixed_loc)

                umi = read_info.seq[umi_fixed_loc: umi_fixed_loc + 16]
                umis_5p.append(umi)
            else:
                umi = read_info.seq[41:57]
                umis_5p.append(umi)
                umi_fixed_locs_5p.append(umi_fixed_loc)
        else:
            putative_bc = read_info.seq[:41]
            putative_bcs_5p.append(putative_bc)
            putative_bc_min_q = _safe_min_q(read_info.q_letter[:41])
            putative_bc_min_qs_5p.append(putative_bc_min_q)

            find_umi_seq = read_info.seq[41:57]
            umi_fixed_loc_re = find_pos(find_umi_seq, FIX6_UMI_5p)

            if umi_fixed_loc_re != -1:
                umi_fixed_loc = umi_fixed_loc_re + 41
                umi_fixed_locs_5p.append(umi_fixed_loc)

                umi = read_info.seq[umi_fixed_loc: umi_fixed_loc + 16]
                umis_5p.append(umi)
            else:
                umi = read_info.seq[41:57]
                umis_5p.append(umi)
                umi_fixed_locs_5p.append(umi_fixed_loc)

    elif BC_fixed_loc > 10:
        # barcode 整体偏后，但通常仍完整
        bc_start = BC_fixed_loc - 10
        bc_end = BC_fixed_loc + 31
        fix5_start = BC_fixed_loc + 16

        putative_bc = read_info.seq[bc_start:bc_end]
        putative_bcs_5p.append(putative_bc)
        putative_bc_min_q = _safe_min_q(read_info.q_letter[bc_start:bc_end])
        putative_bc_min_qs_5p.append(putative_bc_min_q)

        find_umi_seq = read_info.seq[bc_end:bc_end + 16]
        umi_fixed_loc_re = find_pos(find_umi_seq, FIX6_UMI_5p)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + bc_end
            umi_fixed_locs_5p.append(umi_fixed_loc)

            umi = read_info.seq[umi_fixed_loc: umi_fixed_loc + 16]
            umis_5p.append(umi)
        else:
            umi = read_info.seq[bc_end:bc_end + 16]
            umis_5p.append(umi)
            umi_fixed_locs_5p.append(umi_fixed_loc)

    elif BC_fixed_loc < 10 and BC_fixed_loc != -1:
        # barcode 靠左，BC1 可能残缺
        bc_end = BC_fixed_loc + 31

        putative_bc = read_info.seq[:bc_end]
        putative_bcs_5p.append(putative_bc)
        putative_bc_min_q = _safe_min_q(read_info.q_letter[:bc_end])
        putative_bc_min_qs_5p.append(putative_bc_min_q)

        find_umi_seq = read_info.seq[bc_end:bc_end + 16]
        umi_fixed_loc_re = find_pos(find_umi_seq, FIX6_UMI_5p)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re + bc_end
            umi_fixed_locs_5p.append(umi_fixed_loc)

            umi = read_info.seq[umi_fixed_loc: umi_fixed_loc + 16]
            umis_5p.append(umi)
        else:
            umi = read_info.seq[bc_end:bc_end + 16]
            umis_5p.append(umi)
            umi_fixed_locs_5p.append(umi_fixed_loc)

    else:
        # 找不到 FIX6_A_5p
        putative_bcs_5p.append(read_info.seq[:41])
        putative_bc_min_qs_5p.append(_safe_min_q(read_info.q_letter[:41]))

        find_umi_seq = read_info.seq[:80]
        umi_fixed_loc_re = find_pos(find_umi_seq, FIX6_UMI_5p)

        if umi_fixed_loc_re != -1:
            umi_fixed_loc = umi_fixed_loc_re
            umi_fixed_locs_5p.append(umi_fixed_loc)

            umi = read_info.seq[umi_fixed_loc: umi_fixed_loc + 16]
            umis_5p.append(umi)
        else:
            umi_fixed_locs_5p.append(umi_fixed_loc)
            umis_5p.append(umi)

    return (
        read_ids_5p,
        putative_bcs_5p,
        bc_fixed_locs_5p,
        putative_bc_min_qs_5p,
        umis_5p,
        umi_fixed_locs_5p,
    )

####
def norm_seq(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    return "" if s in ("", "NA", "NAN", "NONE") else s

_rc_map = str.maketrans("ACGTN", "TGCAN")
def revcomp(seq):
    seq = seq.upper()
    return seq.translate(_rc_map)[::-1]

def strip_fixed_3p(seq):
    """
    BC3: BC1(10) + AGATC(5) + BC2(10) + GCTACC(6) + BC3(10) = 41
    输出：BC1 + BC2 + BC3 = 30
    """
    s = norm_seq(seq)
    if not s:
        return ""
    if len(s) == 41 and s[10:15] == FIX5_B and s[25:31] == FIX6_A:
        bc1 = s[0:10]
        bc2 = s[15:25]
        bc3 = s[31:41]
        return bc1 + bc2 + bc3
    return ""

def strip_fixed_5p(seq):
    """
    BC5: BC1(10) + CCTTCC(6) + BC2(10) + TGCTG(5) + BC3(10) = 41
    输出：BC1 + BC2 + BC3 = 30
    """
    s = norm_seq(seq)
    if not s:
        return ""
    if len(s) == 41 and s[10:16] == FIX6_A_5p and s[26:31] == FIX5_B_5p:
        bc1 = s[0:10]
        bc2 = s[16:26]
        bc3 = s[31:41]
        return bc1 + bc2 + bc3
    return ""

def norm_bc(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if s in ("", "NA", "nan", "None"):
        return ""
    return s

def build_graph_from_pair_counts(pair_counts_kept, col5="BC5n", col3="BC3n", wcol="support_reads"):
    df = pair_counts_kept[[col5, col3, wcol]].copy()
    df[col5] = df[col5].map(norm_bc)
    df[col3] = df[col3].map(norm_bc)
    df = df[(df[col5] != "") & (df[col3] != "")]
    df[wcol] = df[wcol].astype(int)

    # 无向化 + 聚合同一无向边
    a = df[[col5, col3]].min(axis=1)
    b = df[[col5, col3]].max(axis=1)
    uv = pd.DataFrame({"a": a, "b": b, "w": df[wcol].values})
    edge_agg = uv.groupby(["a", "b"], as_index=False)["w"].sum()

    G = nx.Graph()
    for a, b, w in edge_agg.itertuples(index=False):
        G.add_edge(a, b, weight=int(w))
    return G, edge_agg

def component_category(H: nx.Graph):
    n = H.number_of_nodes()
    m = H.number_of_edges()
    m_self = sum(1 for u,v in H.edges() if u == v)
    m_no_self = m - m_self

    if n == 1 and m == 1 and m_self == 1:
        return "self_only"
    if n == 2 and m_no_self == 1:
        return "pair_only"
    if n == 3 and m_no_self == 3:
        return "triangle_only"
    if n == 4 and m_no_self == 6:
        return "clique4_only"
    return "other"

def collect_components_by_type(G: nx.Graph):
    by = {}
    for comp_nodes in nx.connected_components(G):
        H = G.subgraph(comp_nodes).copy()
        cat = component_category(H)
        by.setdefault(cat, []).append(H)
    return by

# ========= 1) 用 pair_counts_kept 构建 core cells =========
def build_core_cells(pair_counts_kept,
                     allowed_types=("self_only","pair_only","triangle_only","clique4_only"),
                     col5="BC5n", col3="BC3n", wcol="support_reads"):
    G, edge_agg = build_graph_from_pair_counts(pair_counts_kept, col5=col5, col3=col3, wcol=wcol)
    by_type = collect_components_by_type(G)

    core_components = []
    for t in allowed_types:
        core_components.extend(by_type.get(t, []))

    # 给每个 core component 一个 cell_id
    cell_records = []
    barcode2cell = {}
    for i, H in enumerate(core_components, 1):
        cell_id = f"cell_{i:06d}"
        bcs = sorted(H.nodes())
        cell_records.append({
            "cell_id": cell_id,
            "type": component_category(H),
            "n_barcodes": len(bcs),
            "barcodes": bcs
        })
        for bc in bcs:
            barcode2cell[bc] = cell_id

    core_cells_df = pd.DataFrame(cell_records)
    return G, edge_agg, core_cells_df, barcode2cell

# ========= 2) 计算 Top1 + dominance（基于 pair_counts_kept） =========
def compute_top1_dominance(pair_counts_kept, col5="BC5n", col3="BC3n", wcol="support_reads"):
    pc = pair_counts_kept[[col5, col3, wcol]].copy()
    pc[col5] = pc[col5].map(norm_bc)
    pc[col3] = pc[col3].map(norm_bc)
    pc = pc[(pc[col5] != "") & (pc[col3] != "")]
    pc[wcol] = pc[wcol].astype(int)

    # 无向视角展开：barcode -> partner
    long = pd.concat([
        pc.rename(columns={col5:"barcode", col3:"partner", wcol:"w"}),
        pc.rename(columns={col3:"barcode", col5:"partner", wcol:"w"}),
    ], ignore_index=True)

    # 对每个 barcode，按 partner 聚合权重
    bp = long.groupby(["barcode","partner"], as_index=False)["w"].sum()

    # 每个 barcode：sum_all / top1 partner / top1 weight / dominance
    bp_sorted = bp.sort_values(["barcode","w"], ascending=[True, False])
    top1 = bp_sorted.groupby("barcode").head(1).rename(columns={"partner":"top1_partner","w":"top1_w"})
    sum_all = bp.groupby("barcode")["w"].sum().reset_index(name="sum_all")

    dom = sum_all.merge(top1[["barcode","top1_partner","top1_w"]], on="barcode", how="left")
    dom["dominance"] = dom["top1_w"] / dom["sum_all"]

    return dom  # columns: barcode, sum_all, top1_partner, top1_w, dominance

# ========= 3) 将 reads 归到 core cells（规则A + 规则B） =========
def assign_reads_to_cells(df_reads,
                          barcode2cell,
                          dom_table,
                          dominance_min=0.80,
                          bc5_col="BC5n",
                          bc3_col="BC3n",
                          read_id_col="read_id"):
    df = df_reads.copy()
    df["_b5"] = df[bc5_col].map(norm_bc)
    df["_b3"] = df[bc3_col].map(norm_bc)

    # 先判 read 类型
    df["has5"] = df["_b5"] != ""
    df["has3"] = df["_b3"] != ""
    df["read_kind"] = (
        df["has5"].astype(int).astype(str) + "+" + df["has3"].astype(int).astype(str)
    )  # "1+1" paired, "1+0" only5, "0+1" only3, "0+0" none

    # 规则A：直接命中 core barcodes
    df["cell_A_5"] = df["_b5"].map(barcode2cell)
    df["cell_A_3"] = df["_b3"].map(barcode2cell)

    # paired read：如果两端都能命中且不一致，标记冲突（你可以后续决定怎么处理）
    df["cell_A"] = df["cell_A_5"].combine_first(df["cell_A_3"])
    df["cell_A_conflict"] = (df["cell_A_5"].notna()) & (df["cell_A_3"].notna()) & (df["cell_A_5"] != df["cell_A_3"])

    # 规则B：对 single-end（且规则A没命中）尝试用 top1+dominance 吸附
    dom = dom_table.set_index("barcode")

    def try_absorb(bc):
        if bc == "" or bc not in dom.index:
            return None
        row = dom.loc[bc]
        if pd.isna(row["top1_partner"]):
            return None
        if float(row["dominance"]) < dominance_min:
            return None
        # top1 partner 必须在 core cells
        return barcode2cell.get(row["top1_partner"], None)

    need_B = df["cell_A"].isna() & (df["has5"] ^ df["has3"])  # 单端且A没命中
    df.loc[need_B & df["has5"], "cell_B"] = df.loc[need_B & df["has5"], "_b5"].map(try_absorb)
    df.loc[need_B & df["has3"], "cell_B"] = df.loc[need_B & df["has3"], "_b3"].map(try_absorb)

    # 最终 cell_id：优先 A，其次 B
    df["cell_id"] = df["cell_A"].combine_first(df["cell_B"])

    # 一些统计
    stats = {
        "n_total_reads": len(df),
        "n_paired": int((df["read_kind"]=="1+1").sum()),
        "n_single_5": int((df["read_kind"]=="1+0").sum()),
        "n_single_3": int((df["read_kind"]=="0+1").sum()),
        "n_assigned_A": int(df["cell_A"].notna().sum()),
        "n_assigned_B_only": int((df["cell_A"].isna() & df["cell_B"].notna()).sum()),
        "n_unassigned": int(df["cell_id"].isna().sum()),
        "n_conflict_paired_A": int(df["cell_A_conflict"].sum()),
    }

    # 清理临时列（你也可以保留用于debug）
    df = df.drop(columns=["_b5","_b3"], errors="ignore")

    return df, stats

def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def iter_fastq(handle):
    """yield (header, seq, plus, qual) each with trailing \n"""
    while True:
        h = handle.readline()
        if not h:
            return
        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        if not q:
            return
        yield h, s, p, q

def extract_reads_and_filter_df_by_raw(
    df: pd.DataFrame,
    raw_fastq_gz: str,
    out_fastq_gz: str,
    read_id_col: str = "read_id",
    remove_found: bool = True,
    return_missing_ids: bool = False,   # True 时会把 missing_ids 也返回（可能很大）
):
    """
    1) 从 raw_fastq_gz 中提取 df[read_id_col] 对应的 reads，写入 out_fastq_gz（gzip）
    2) df 中如果某 read_id 没在 raw 里找到：从 df 删除该行
    返回:
      df_kept, stats
      其中 stats 包含 dropped_rows 等
    """
    if read_id_col not in df.columns:
        raise KeyError(f"df 中没有列: {read_id_col}")

    # 标准化 df 里的 read_id
    df2 = df.copy()
    df2[read_id_col] = df2[read_id_col].astype(str).str.strip()
    df2 = df2[df2[read_id_col].notna() & (df2[read_id_col] != "")].copy()

    # 目标集合（唯一）
    id_set = set(df2[read_id_col].tolist())
    n_target_unique = len(id_set)
    n_target_rows = len(df2)

    found_ids = set()
    found = 0
    scanned = 0

    with open_maybe_gz(raw_fastq_gz, "rt") as fin, gzip.open(out_fastq_gz, "wt") as fout:
        for h, s, p, q in iter_fastq(fin):
            scanned += 1
            rid = h[1:].strip().split()[0]  # '@'后到空格前

            if rid in id_set:
                fout.write(h); fout.write(s); fout.write(p); fout.write(q)
                found += 1
                found_ids.add(rid)

                if remove_found:
                    id_set.remove(rid)
                    if not id_set:
                        break

    # 用 found_ids 过滤 df（删掉没找到的行）
    df_kept = df2[df2[read_id_col].isin(found_ids)].copy()

    # 统计：注意 df2 可能有重复 read_id，因此 dropped_rows 用行数算
    dropped_rows = int(n_target_rows - len(df_kept))
    missing_unique = int(n_target_unique - len(found_ids))

    stats = {
        "target_rows": int(n_target_rows),           # df 中参与匹配的行数（含重复 read_id）
        "target_unique_ids": int(n_target_unique),   # 唯一 read_id 数
        "found_unique_ids": int(len(found_ids)),     # 找到的唯一 read_id 数
        "written_reads": int(found),                 # 写入 fastq 的 reads 数（remove_found=True 时一般=found_unique_ids）
        "scanned_reads": int(scanned),
        "dropped_rows": int(dropped_rows),           # df 被删掉的行数
        "missing_unique_ids": int(missing_unique),   # 没找到的唯一 read_id 数
        "out_fastq_gz": out_fastq_gz,
    }

    if return_missing_ids:
        stats["missing_ids"] = sorted(list(set(df2[read_id_col]) - found_ids))

    return df_kept, stats

def is_missing(x):
    # True 表示“空/缺失”
    if pd.isna(x):
        return True
    s = str(x).strip()
    return s == "" or s.upper() in {"NA", "NAN", "NONE"}

def umi_missing_stats(df, col3="putative_umi", col5="putative_umi_5p"):
    m3 = df[col3].map(is_missing)  # True=missing
    m5 = df[col5].map(is_missing)

    n = len(df)

    # 你要的三类 + 额外给出“分别不为空”
    both_present = (~m3 & ~m5).sum()
    only_3p = (~m3 &  m5).sum()
    only_5p = ( m3 & ~m5).sum()
    both_missing = ( m3 &  m5).sum()

    out = pd.DataFrame({
        "category": [
            "putative_umi present",
            "putative_umi_5p present",
            "both present",
            "only putative_umi present",
            "only putative_umi_5p present",
            "both missing",
        ],
        "count": [
            (~m3).sum(),
            (~m5).sum(),
            both_present,
            only_3p,
            only_5p,
            both_missing,
        ],
    })
    out["ratio"] = out["count"] / n
    return out


