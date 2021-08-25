#!/usr/bin/env python3

# Consume results from a dataset test, building a table for analysis and making plots

import os
import subprocess
import shutil
import argparse

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from itertools import takewhile
from matplotlib import rc
import os
import sys
import signal
import random
from math import log

def signal_handler(signal, frame):
    print("\nprogram exiting gracefully")
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

# Style options
sns.set()

sns.set(font="Linux Biolinum O", font_scale=1.5)
# rc('text', usetex=True)

keenan_purple = "#1b1f8a"
alt_color = "#1f8a1b"
sad_color = "#ff0000"

def ensure_dir_exists(d):
    if not os.path.exists(d):
        os.makedirs(d)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('dir', help="directory with result")
    parser.add_argument('--record_files', action='store_true')
    parser.add_argument('--max_files', type=int, default = -1, help="read in at most max_files many of the input files")
    parser.add_argument('--name', default = "", help="data set name to use in captions")
    parser.add_argument('--merged_files', default = "", help="csv of dataset")
    parser.add_argument('--plot_runtime', action='store_true', help="plot parameterization runtime")
    args = parser.parse_args()

    plot_time_vs_mesh_size    = args.plot_runtime

    # Place results here
    out_dir = os.path.join(args.dir, "analysis")
    ensure_dir_exists(out_dir)

    # Parse in the individual CSVs
    frames = []
    bad_frames = []
    bad_files = []

    synthetic_timeout = 2000

    if args.merged_files:
        df = pd.read_csv(args.merged_files, header=0, sep=',',engine='python')
        # drop rows whose times are too long
        bad_rows = df[df['duration'] > synthetic_timeout ].index
        df.drop(bad_rows, inplace=True)
        bad_df = None
        print(df)
        print("Done loading df")
    else:
        data_files = os.listdir(args.dir)
        if (args.max_files > 0):
            data_files = os.listdir(args.dir)[:args.max_files]
        for i,f in enumerate(data_files):
            if not f.endswith(".tsv"):
                continue

            if i % 10 == 0:
                # print("loading " + str(f))
                wdth = 50
                progress = int(i / len(data_files) * wdth)
                bar = "[" + (progress  * "=") + ((wdth-progress)*" ") + "]"
                print(f"\rLoading Data {bar} ({int(progress / wdth * 100)}%)", end="", flush=True)
                # print("loading " + str(f))
            try :
                frame = pd.read_csv(os.path.join(args.dir,f), header=0, sep='\t',engine='python')

                # remove leading/trailing whitespace from column names
                frame = frame.rename(columns=lambda x:x.strip())

                # replace inf w/ np.inf, allowing for spaces around "inf"
                frame = frame.replace(to_replace=r"^\s*inf\s*$", value=np.inf, regex=True)

                # replace nan w/ np.nan, allowing for spaces around "nan"
                frame = frame.replace(to_replace=r"^\s*nan\s*$", value=np.nan, regex=True)

                # replace " false" w/ False, allowing for spaces around "false"
                frame = frame.replace(to_replace=r"^\s*false\s*$", value=False, regex=True)
            except:
                frame = pd.DataFrame()
            if frame.empty or "name" not in frame.columns:
                bad_files.append(f)
            elif not isinstance(frame.at[0, "nVertices"], (int,float, np.integer)) or not isinstance(frame.at[0, "duration"], (int, float, np.integer)):
                print(f"Vertex type {frame.at[0, 'nVertices']} : {type(frame.at[0, 'nVertices'])}")
                print(f"duration type {frame.at[0, 'duration']} : {type(frame.at[0, 'duration'])}")
                bad_files.append(f)
            # elif frame.shape[1] != 20 or pd.isnull(frame.at[0, 'nvertices']):
            elif pd.isnull(frame.at[0, 'nVertices']):
                bad_frames.append(frame.head(1))
            else:
                frames.append(frame)
        print("\n")

        # Make a mega-frame
        if frames:
            df = pd.concat(frames, ignore_index=True)
        else:
            df = pd.DataFrame()

        # Save it
        df.to_csv(os.path.join(out_dir, "merged_results.csv"), index=False)

        bad_df = None
        if(len(bad_frames) > 0):
            bad_df = pd.concat(bad_frames, ignore_index=True)

        # drop rows whose times are too long
        bad_rows = df[df['duration'] > synthetic_timeout].index
        df.drop(bad_rows, inplace=True)

        if bad_files:
            print(f"\n\n === bad files({len(bad_files)}) === ")
            print(bad_files)
        else:
            print("\n\n === all files parsed :) ")

        if bad_df is None:
            print("\n\n === no incomplete records :) ")
            filtered_bad_names = []
        else:
            print("\n\n === incomplete records === ")
            print(bad_df)

            bad_names = bad_df["name"].to_list()

            filtered_bad_names = [int(c) if len(c) else None for c in (''.join(takewhile(str.isdigit, str(x) or "")) for x in bad_names)]

        # print(bad_names)
        # print(filtered_bad_names)
        with open(f"bad_meshes.txt", "w") as f:
            f.write("\n".join(map(str, [x[:-4] for x in bad_files])))
            f.write("\n")
            f.write("\n".join(map(str, filtered_bad_names)))


        if args.record_files:
            with open("good_meshes.txt", "w") as f:
                f.write("\n".join(map(str,df["name"])))

            df['relevant time'] = df['duration'] - df['cone placement time']
            slow_meshes=df[(df["relevant time"] >= 1000)]
            with open("slow_meshes.txt", "w") as f:
                f.write("\n".join(map(str,slow_meshes["name"])))

            noninjective_meshes=df[df["nFlippedTriangles"] + df["nZeroAreaTriangles"] > 0]
            with open("noninjective_meshes.txt", "w") as f:
                f.write("\n".join(map(str,noninjective_meshes["name"])))


    """
    # Filtering, etc
    df = df[df["minCornerAngle"] > 1e-6]

    # Compute some extra values
    df["flip_timeFrac"] = df["flip_flipTime"] / df["flip_totalTime"]

    # Print the table to the terminal
    """
    print(df)
    print(df.columns)

    n_bad_df = 0 if bad_df is None else bad_df.shape[0]

    n_total = len(df) + n_bad_df + len(bad_files)
    if args.name =="Thingi10k":
        n_total = 32744 # number of total thingi10k meshes
    else:
        print(f"args.name: {args.name}")

    print(f"\tfailed to terminate successfully within the time limit on {n_bad_df+len(bad_files)} meshes")

    n = df.shape[0]
    df['injective'] = df.apply(lambda x : (x['nFlippedTriangles'] == 0) and (x['nZeroAreaTriangles']==0), axis=1)
    print(f"Achieved local injectivity on {df['injective'].sum()} of {n} meshes ({df['injective'].sum()/n*100.}% of successes, {df['injective'].sum()/n_total*100.}% of all meshes)")

    print(df[df['injective']==False])

    mean_runtime = df['duration'].mean()
    print(f"mean runtime: {mean_runtime}")

    df_not_tiny = df[df['nVertices'] > 1000]
    mean_not_tiny_runtime = df_not_tiny['duration'].mean()
    print(f"mean not tiny runtime: {mean_not_tiny_runtime}")


    ############################################################
    #                   Plots
    ############################################################

    name = args.name

    ### Plot Time vs Mesh Size
    if plot_time_vs_mesh_size:
        scale = 2
        fig, ax = plt.subplots(figsize=(3.377 * scale,1.75 * scale))

        ax.scatter(df["nVertices"], df["duration"], color=keenan_purple, edgecolors='none',alpha=0.5)

        ax.set_xlabel("number of mesh vertices")
        ax.set_ylabel("duration (s)")
        ax.set_title(f"{name} Time vs Mesh Size")

        fig.savefig(os.path.join(out_dir,'time.pdf'))

        plt.show()
        plt.close()

if __name__ == "__main__":
    main()
