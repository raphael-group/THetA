from FileIO import *

filename,results,n,k,tau,directory,prefix,max_normal,bound_heuristic,normal_bound_heuristic, heuristic_lb, heuristic_ub, num_processes, bounds_only, multi_event, force, get_values, interval_selection, num_intervals, read_depth_file, graph_format = parse_arguments(silent = True)

#runN3 = os.path.join(directory, prefix+".RunN3.bash")
prefix = os.path.join(directory,prefix)

print prefix

