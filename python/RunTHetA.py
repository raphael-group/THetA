 ###
 # Copyright 2012, 2013, 2014, 2015, 2017 Brown University, Providence, RI.
 #
 # All Rights Reserved
 # 
 # Permission to use this software, and any documentation, for non-commercial academic research
 # purposes only is hereby granted with the following terms and conditions:
 # 
 # (1) the above copyright notice and this permission notice shall be preserved in all instances
 # of the software and in any supporting documentation;
 # 
 # (2) the name of Brown University shall not be used in advertising or publicity pertaining 
 # to the use of the software without specific, written prior permission;
 # 
 # (3) the rights granted herein are individual and personal to the recipient and may not be
 # sublicensed or distributed to any third party without specific, written prior permission; and
 # 
 # (4) the permitted user acknowledges that all commercial rights are licensed to Medley
 # Genomics, Inc., and any inquiries related to commercial use shall be directed to Medley
 # Genomics, Inc.
 # 
 # BROWN UNIVERSITY PROVIDES THIS SOFTWARE AND ANY DOCUMENTATION
 # "AS IS" AND DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE
 # AND ANY DOCUMENTATION, INCLUDING ALL IMPLIED WARRANTIES OF
 # MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE. IN NO
 # EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 # CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 # LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 # NEGLIGENCE OR OTHER ACTION BASED ON ANY OTHER LEGAL THEORY,
 # ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 # SOFTWARE.
 #
 # @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael, Gryte Satas and
 # Alex Ashery
 ###

from SelectIntervals import *
from FileIO import *
from DataTools import *
from Misc import *
from Enumerator import Enumerator
from Optimizer import Optimizer
from TimeEstimate import *
from CalcAllC import *
from ModelSelection import *
from plotResults import *
from RunBAFModel import run_BAF_model, calculate_BAF, generate_pi, calculate_interval
from ClusteringBAF import clustering_BAF, group_to_meta_interval

from multiprocessing import JoinableQueue, Queue, Process, Array, current_process
import os

def process_loop(queue, opt, returnQueue, sorted_index, get_values):
    """
    Code that each child process executes. Repeatedly pops of new C
    values from queue until it reaches an exit signal. Then puts its results
    on the return queue and finishes

    Arguments:
    queue (multiprocessing.Queue): Task queue, containing C matrices
    opt (Optimizer): instance of an optimizer
    returnQueue (multiprocessing.Queue): Queue to put results in
    sorted_index (list): Array containing ordering information for sorting
    """
    min_likelihood = float('inf')
    best = []
    if get_values: solns = []
    while True:
        C = queue.get()
        if C is 0:
            returnQueue.put(best)
            break
        soln = opt.solve(C)
        if soln is not None:
            (mu, likelihood,vals) = soln
            if get_values: solns.append((C,mu,likelihood,vals))
            if isClose([likelihood],[min_likelihood]):
                C_new = reverse_sort_C(C,sorted_index)
                vals = reverse_sort_list(vals, sorted_index)
                best.append((C_new, mu, likelihood, vals))
            elif likelihood < min_likelihood:
                C_new = reverse_sort_C(C,sorted_index)
                vals = reverse_sort_list(vals, sorted_index)
                best = [(C_new, mu, likelihood, vals)]
                min_likelihood = likelihood

    if get_values:
        with open(pre+"."+"values" + str(current_process().name),'w') as f:
            for C,mu,likelihood,vals in solns:
                m,n = C.shape
                stringC = "".join((str(int(C[i][1])) for i in range(m)))
                valsStr = " ".join((str(v) for v in vals))
                f.write(stringC+"\t"+str(mu[0])+"\t"+str(likelihood)+"\t"+valsStr+"\n")


def start_processes(max_processes, queue, opt, returnQueue, sorted_index, get_values):
    """
    Starts a max_processes number of processes, and starts them
    """
    processes = [Process(target=process_loop, args=(queue, opt, returnQueue,\
                sorted_index, get_values), name=i+1) for i in range(max_processes-1)]
    for p in processes:
        p.daemon = True
        p.start()
    return processes

def find_mins(best):
    """
    Takes a the list of "best" C,mu pairs returned by each process and finds
    the ones with the minimum likelihood
    """
    min_likelihood = float('inf')
    true_best = []
    for solns in best:
        if len(solns) == 0: continue
        likelihood = solns[0][2]
        if isClose([min_likelihood], [solns[0][2]]):
            true_best += solns
        elif likelihood < min_likelihood:
            min_likelihood = likelihood
            true_best = solns
    return true_best

def do_optimization(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
            max_normal, sorted_index, max_processes, multi_event, get_values):
    """
    Performs the optimization for the given parameters with max_proccesses
    number of processes
    Returns a list of the best C matrices and associated mu values
    and likelihoods
    """
    enum = Enumerator(n, m, k, tau, lower_bounds, upper_bounds, multi_event)
    opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)
    MAX_QUEUE_SIZE = int(10E6)
    try:
        queue = Queue(MAX_QUEUE_SIZE) #Task queue for the processes
    except OSError:
        MAX_QUEUE_SIZE = 2**15-1
        queue = Queue(MAX_QUEUE_SIZE)

    returnQueue = Queue(MAX_QUEUE_SIZE) #Shared queue for processes to return results

    processes = start_processes(max_processes, queue, opt, returnQueue, \
                sorted_index, get_values)

    # fix problem with missing first matrix
    #C = enum.generate_next_C()
    C=enum._C_to_array()
    count = 0
    while C is not False:
        count += 1
        queue.put(C, True)
        C = enum.generate_next_C()
    if count == 0:
        print "Error: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
        sys.exit(1)

    # Send STOP signal to all processes
    for i in range(max_processes-1):
        queue.put(0)

    best = []
    for i in range(len(processes)):
        item = returnQueue.get()
        best.append(item)

    for p in processes:
        p.join()

    best = find_mins(best)
    return best

def do_optimization_single(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
            max_normal, sorted_index, multi_event, get_values):
    """
    Performs the optimization for the given parameters with a single process
    Returns a list of the best C matrices and associated mu values
    and likelihoods
    """

    enum = Enumerator(n, m, k, tau, lower_bounds, upper_bounds, multi_event)
    opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)
    min_likelihood = float("inf")
    best = []
    count = 0

    #fix missing first matrix problem
    C=enum._C_to_array()
    #C = enum.generate_next_C()
    if get_values: solns = []
    while C is not False:
        count += 1
        soln = opt.solve(C)
        if soln is not None:
            (mu, likelihood,vals) = soln

            if get_values: solns.append((C,mu,likelihood))
            if isClose([likelihood],[min_likelihood]):
                C_new = reverse_sort_C(C,sorted_index)
                vals = reverse_sort_list(vals, sorted_index)
                best.append((C_new, mu, likelihood, vals))
            elif likelihood < min_likelihood:
                C_new = reverse_sort_C(C,sorted_index)
                vals = reverse_sort_list(vals, sorted_index)
                best = [(C_new, mu, likelihood, vals)]
                min_likelihood = likelihood

        C = enum.generate_next_C()

    if get_values:
        with open(pre+"."+"likelihoods",'w') as f:
            for C,mu,likelihood in solns:
                m,n = C.shape
                stringC = "".join((str(int(C[i][1])) for i in range(m)))
                f.write(stringC+"\t"+str(mu[0])+"\t"+str(likelihood)+"\n")

    if count == 0:
        print "Error: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
        sys.exit(1)
    return best

def best_near_max_contamination(best, max_normal):
    for C, mu, likelihood, vals in best:
        if abs(max_normal-mu[0]) < .01: return True
    return False

def get_clustering_args(tumorfile, normalfile, filename, num_processes, m, tumorCounts, normCounts):
    tumorData = read_snp_file(tumorfile)
    normalData = read_snp_file(normalfile)
    chrmsToUse, intervalData = read_interval_file_BAF(filename)
    minSNP = 10
    gamma = 0.05
    print "Calculating BAFs"
    tumorBAF, normalBAF, tumorData, normalData = calculate_BAF(tumorData, normalData, chrmsToUse, minSNP, gamma, num_processes)

    pi = generate_pi(intervalData)
    SNPToIntervalMap = [calculate_interval(pi, snp[0], snp[1]) for snp in tumorData]
    meanBAFs = [0 for i in range(m)]

    numSNPs = [0 for i in range(m)]
    for i in range(len(SNPToIntervalMap)):
        mapping = SNPToIntervalMap[i]
        if mapping is None: continue
        meanBAFs[mapping] += abs(tumorBAF[i] - 0.5)
        numSNPs[mapping] += 1.0
    meanBAFs = map(lambda (num, denom): num / denom if denom > 0 else -1, zip(meanBAFs, numSNPs))

    corrRatio = []
    tTotal = float(sum(tumorCounts))
    nTotal = float(sum(normCounts))
    for i in range(m):
        tCount = float(tumorCounts[i])
        nCount = float(normCounts[i])
        if nCount == 0 or meanBAFs[i] == -1:
            corrRatio.append(-1)
            meanBAFs[i] = -1
        else:
            corrRatio.append((tCount / tTotal) / (nCount / nTotal))

    chrms, starts, ends = zip(*intervalData)
    intervals = zip(chrms, starts, ends, tumorCounts, normCounts, corrRatio, meanBAFs, numSNPs)

    intervalsByChrm = [[] for i in range(24)]
    missingData = []
    for i, interval in enumerate(intervals):
        if interval[5] == -1 or interval[6] == -1:
            interval = list(interval)
            interval.append(i)
            missingData.append(interval)
        else:
            chrm = interval[0]
            intervalsByChrm[chrm].append(interval)

    intervals = intervalsByChrm

    return intervals, missingData, corrRatio, meanBAFs, tumorData, normalData, tumorBAF, normalBAF, chrmsToUse, intervalData

def main():
    ###
    #  Read in arguments and data file
    ##
    args = parse_arguments()
    print "Reading in query file..."
    intervals = read_interval_file(args[0])

    if args[2] != None:
        # Run for a specific n. If n > 2, args[1] should contain the results file
        # Backwards compatibility
        resultsfile, boundsfile = run_fixed_N(args[2], args, intervals, args[1])
    else:
        resultsfile2, boundsfile2 = run_fixed_N(2, args, intervals)
        intervals = read_interval_file(boundsfile2)
        resultsfile3, boundsfile3 = run_fixed_N(3, args, intervals, resultsfile2)

        ModelSelection(args[0],resultsfile2, resultsfile3)


def run_fixed_N(n, args, intervals, resultsfile=None):
    (filename, results, N, k, tau, directory, prefix, max_normal, bound_heuristic, \
        normal_bound_heuristic,heuristic_lb, heuristic_ub, num_processes, \
        bounds_only, multi_event, force, get_values, choose_intervals, num_intervals, \
        read_depth_file, graph_format, runBAF, ratio_dev, min_frac,\
        tumorfile, normalfile, noClustering) = args

    lengths, tumorCounts, normCounts, m, upper_bounds, lower_bounds = intervals

    global pre
    pre = prefix

    ###
    #   Determine is sample has enough copy number aberrations to run
    ###
    frac = determine_frac_copy_num(normCounts, tumorCounts, lengths, ratio_dev)
    print "Frac with potential copy numbers:", frac
    if frac < min_frac:
        print "ERROR: This sample does not have enough large copy number aberrations to be a good candidate for tumor composition estimation using THetA.  See --RATIO_DEVIATION and --MIN_FRAC flags to modify how the potential presence of large copy number aberrations is determined.  Exiting..."
        exit(1)

    doClustering = tumorfile is not None and normalfile is not None and not noClustering

    ###
    # Setup if we will do clustering of intervals prior to running
    ###
    if doClustering:
        intervals, missingData, corrRatio, meanBAFs, tumorData, normalData, tumorBAF, normalBAF, chrmsToUse, intervalData = get_clustering_args(tumorfile, normalfile, filename, num_processes, m, tumorCounts, normCounts)

        #original clustering code
        lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters, clusterMeans, normalInd = clustering_BAF(n, intervals=intervals, missingData=missingData, prefix=prefix, outdir=directory, numProcesses=num_processes)

        origM, origLengths, origTumor, origNormal, origUpper, origLower = (m, lengths, tumorCounts, normCounts, upper_bounds, lower_bounds)

        intervalMap, lengths, tumorCounts, normCounts, lower_bounds, upper_bounds = \
        group_to_meta_interval(lengths, tumorCounts, normCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters)

        m = len(lengths)

        cluster_scores = score_clusters(intervalMap, origLengths, corrRatio, meanBAFs, m)


    ###
    #   Automatically Select Intervals
    #   note: This is the default behavior
    ###
    if choose_intervals:

        if doClustering:
            print "Selecting meta-intervals..."
            allM, allLengths, allTumor, allNormal, allUpperBounds, allLowerBounds = (origM, origLengths, origTumor, origNormal, origUpper, origLower)
        else:
            print "Selecting intervals..."
            allM, allLengths, allTumor, allNormal, allUpperBounds, allLowerBounds = (m, lengths, tumorCounts, normCounts, upper_bounds, lower_bounds)

        if n == 2:

            if doClustering:

                order, lengths, tumorCounts, normCounts, lower_bounds, upper_bounds = \
                select_meta_intervals_n2(lengths, tumorCounts, normCounts, m, k, force, num_intervals, cluster_scores, lower_bounds, upper_bounds)

            else:
                # 07-09-15 - Make so we can set specific bounds
                if lower_bounds is None or upper_bounds is None:
                    order, lengths, tumorCounts, normCounts = select_intervals_n2(lengths, tumorCounts, normCounts, m, k, force, num_intervals)
                    upper_bounds = None
                    lower_bounds = None
                else:
                    order, lengths, tumorCounts, normCounts, lower_bounds, upper_bounds = \
                    select_intervals_n2(lengths, tumorCounts, normCounts, m, k, force, num_intervals, lower_bounds, upper_bounds)

        elif n == 3:

            if doClustering:

                order, lengths, tumorCounts, normCounts, lower_bounds, upper_bounds = \
                select_meta_intervals_n3(lengths, tumorCounts, normCounts, m, k, force, num_intervals, cluster_scores, lower_bounds, upper_bounds)

            elif resultsfile is None:
                print "ERROR: No results file supplied. Unable to automatically select intervals for n=3 without results of n=2 analysis. See --RESULTS flag, or --NO_INTERVAL_SELECTION to disable interval selection. Exiting..."
                exit(1)
            else:
                # Need to read in original file, bounds file and results file. Original file needed because copy numbers are based on
                copy = read_results_file(resultsfile)
                order, lengths, tumorCounts, normCounts, upper_bounds, lower_bounds, copy = select_intervals_n3(lengths, tumorCounts, normCounts, m, upper_bounds, lower_bounds, copy, tau, force, num_intervals)

        m = len(order)


    sum_r = sum(tumorCounts)
    sum_rN = sum(normCounts)
    set_total_read_counts(sum_r, sum_rN)
    ###
    #  Process/sort read depth vectors and calculate bounds if necessary
    ###
    print "Preprocessing data..."


    r,rN,sorted_index = sort_r(normCounts,tumorCounts)



    if normal_bound_heuristic is not False:
        upper_bounds,lower_bounds = calculate_bounds_normal_heuristic( \
            normal_bound_heuristic, heuristic_lb, heuristic_ub, r, rN, m, k)
    elif bound_heuristic is not False or upper_bounds is None and lower_bounds is None:
        if bound_heuristic is False: bound_heuristic = 0.5
        upper_bounds,lower_bounds = calculate_bounds_heuristic(float(bound_heuristic),\
             r, rN, m, tau, k)
    else:
        if upper_bounds is not None: upper_bounds = sort_by_sorted_index(upper_bounds,\
            sorted_index)
        if lower_bounds is not None: lower_bounds = sort_by_sorted_index(lower_bounds,\
            sorted_index)



    ###Bounds files in their original orders
    ub_out = reverse_sort_list(upper_bounds, sorted_index)
    lb_out = reverse_sort_list(lower_bounds, sorted_index)

    #Need to un-meta cluster before writing bounds file
    if doClustering:
        ub_out, _ = un_meta_cluster_bounds(ub_out, order, intervalMap)
        meta_order = order
        lb_out, order = un_meta_cluster_bounds(lb_out, order, intervalMap)


    if choose_intervals:
        boundsfile = write_out_bounds(directory, prefix, filename, ub_out, lb_out, n, order)
    else:
        boundsfile = write_out_bounds(directory, prefix, filename, ub_out, lb_out, n)

    if bounds_only: sys.exit(0)



    enum = time_estimate(n,m,k,tau,lower_bounds,upper_bounds,r,rN,max_normal,sorted_index, num_processes, multi_event, force)
    ###
    #  Initialize optimizer and enumerator
    ###
    print "Performing optimization..."

    if num_processes == 1:
        best = do_optimization_single(n,m,k,tau,lower_bounds,upper_bounds,
            r,rN,max_normal,sorted_index, multi_event, get_values)
    else:
        best = do_optimization(n, m, k, tau, lower_bounds, upper_bounds, r, rN,\
                max_normal, sorted_index, num_processes, multi_event, get_values)
    if best == []:
        print "ERROR: Maximum Likelihood Solution not found within given bounds."
        exit(1)

    if n == 2 and best_near_max_contamination(best, max_normal):
        print "WARNING: At least one of the top solutions is near the upper bound on normal contamination. Further analysis may required as the sample likely falls into one of the following categories:\n\t1. This sample has high normal contamination. Consider re-running with an increased normal contamination upper bound. See --MAX_NORMAL option\n\t2. This sample may not satisfy the assumption that most of the tumor genome retains the normal expected copynumber (e.g. a genome duplication event has occurred). See THetA optional parameters in changing the expected copy number.\n\t3. This sample may not be a good candidate for THetA analysis (i.e. does not contain large copy number aberrations that distinguish populations)."
    r = reverse_sort_list(r, sorted_index)
    rN = reverse_sort_list(rN, sorted_index)


    if doClustering:
        if n == 2:
            best, r, rN = un_meta_cluster_results_N2(best, meta_order, intervalMap, allTumor, allNormal)
        elif n == 3:
             best, r, rN = un_meta_cluster_results_N3(best, meta_order, intervalMap, allTumor, allNormal, n)

        #r = XXX
        #rN = XXX

    if choose_intervals:
        if n == 2:
            best = calc_all_c_2(best, r, rN, allTumor, allNormal, order)
        elif n == 3 and not multi_event:
            best = calc_all_c_3(best, r, rN, allTumor, allNormal, order)
        else:
            best = calc_all_c_3_multi_event(best, r, rN, allTumor, allNormal, order)

        #lko 7-6-15 Fix multiple solutions to only return the solution with the overal min NLL
        best = find_mins(best)


    #run BAF model on results to determine most likely solution
    if runBAF and tumorfile is not None and normalfile is not None:
        if len(best) != 1:
            BAFprefix = prefix # + ".preliminary"
            resultsfile = write_out_result(directory, BAFprefix, best, n)

            resultsFile = BAFprefix + ".n"+str(n)+".results"
            resultsPath = os.path.join(directory, resultsFile)
            try:
                run_BAF_model(resultsPath, tumor=tumorData, normal=normalData, normalBAF=normalBAF, tumorBAF=tumorBAF, chrmsToUse=chrmsToUse,intervals=intervalData, prefix=prefix + ".n" + str(n), directory=directory, numProcesses=num_processes)
                #run_BAF_model(resultsPath, prefix=prefix + ".n" + str(n), directory=directory, numProcesses=num_processes)
            except IOError:
                print "ERROR: Invalid locations for tumor and normal SNP files. The BAF model will not be run. You can try running the BAF model again directly from the runBAFModel.py script."
        else:
            resultsfile = write_out_result(directory, prefix, best, n)
    elif runBAF and (tumorfile is None or normalfile is None):
        print "ERROR: Need file location for tumor and normal SNP files to run the BAF model. The BAF model will not be run. You can try running the BAF model again directly from the runBAFModel.py script."
        resultsfile = write_out_result(directory, prefix, best, n)
    else:
        resultsfile = write_out_result(directory, prefix, best, n)

    ###
    # Make Results Plots
    ###
    print "Plotting results as a " + graph_format + "..."
    #plot_results(directory, filename, prefix, n)
    plot_results(directory, filename, prefix, read_depth_file, n, graph_format)

    if n == 2: write_out_N3_script(directory, prefix, filename)

    return resultsfile, boundsfile

import time

if __name__ == '__main__':
    main()
