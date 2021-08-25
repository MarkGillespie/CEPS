## Running a benchmark
`run_benchmark.py` runs the cone flattening code on a dataset that you provide and records the algorithm's performance.
The script produces a `tsv` file recording some performance statistics for each mesh in the input. These files can then be analyed by `summarize_results.py`.

|flag | purpose|
| ------------- |-------------| 
|`--dataset_dir=/path/to/meshes`| Directory of mesh files to run on (required). |
|`--output_dir=/path/to/output/dir`| Directory to store results (required). |
|`--good_list=meshes_to_use.txt`| File of meshes which should be used. By default, all meshes are used. |
|`--bad_list=meshes_to_skip.txt`| File of meshes which should be excluded. |
|`--n_threads=1` | Number of threads to run on (default=1). |
|`--timeout=600` | Timeout in seconds (default=600). |
|`--max_meshes=1000` | Maximum number of meshes to process. By default, all meshes are used. |
|`--use_ffield_cones` | Use cones from an MPZ-style `.ffield` file (must have the same base name as the mesh file, and be in the same directory. |
|`--save_parameterized_meshes` | Save parameterizations as well as performance statistics. |

## Summarizing benchmark results
Run `process_results.py your_output_dir` to summarize the results and plot the algorithm runtime (saved to `your_output_dir/analysis`) along with a `csv` containing the statistics from all meshes.

|flag | purpose|
| ------------- |-------------| 
|`--name=DatasetName` | Dataset name to use in figure captions. |
|`--plot_runtime` | Plot the algorithm runtime. |
|`--merged_files` | Name of a `csv` to read, instead of reading in individual mesh records. |
