# Helper class to run tasks in multiple processes

import os
import subprocess
import shutil
import random
import argparse
import sys

from threading import Thread
from queue import Queue


# Simple task queue
class ShellTaskQueue(Queue):

    def __init__(self, nWorkers=1, timeout=99999999999):
        Queue.__init__(self)
        self.nWorkers = nWorkers
        self.timeout = timeout
        self.task_count = 0

        for i in range(self.nWorkers):
            t = Thread(target=self.worker)
            t.daemon = True
            t.start()

    def add_task(self, cmd_str):
        self.put((self.task_count, cmd_str))
        self.task_count += 1

    def worker(self):
        while True:
            id, cmd_str = self.get()
            print("running task {} / {}    {}\n".format(id+1, self.task_count, cmd_str))

            with subprocess.Popen(cmd_str, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'), preexec_fn=os.setsid) as process:

                try:
                    output = process.communicate(timeout=self.timeout)[0]
                except subprocess.TimeoutExpired:
                    print("  Timeout on {} :(".format(cmd_str))
                    os.killpg(process.pid, subprocess.signal.SIGINT)  # send signal to the process group
                    output = process.communicate()[0]
                except Exception as e:
                    print("  EXCEPTION ON {} !!!!".format(cmd_str))
                    print(str(e))


            self.task_done()

# location of flip binaries, assumed relative to this file
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "bin"))

def ensure_dir_exists(d):
    if not os.path.exists(d):
        os.makedirs(d)

def parse_file_list(f):
    names = set()
    for line in open(f, 'r').readlines():
        if line[0] == '#': continue

        raw_name = line.strip()
        name , _ = os.path.splitext(os.path.basename(raw_name))
        names.add(name)

    return names

def main():
    parser = argparse.ArgumentParser()

    # Build arguments
    parser.add_argument('--dataset_dir', type=str, required=True, help='directory of mesh files to run on')
    parser.add_argument('--output_dir', type=str, required=True, help='where to put results')

    parser.add_argument('--good_list', type=str, help='file with meshes which should be used (optional, otherwise all used. ignores paths and extensions.)')
    parser.add_argument('--bad_list', type=str, help='file with meshes which should not be used, even if on good list (optional, otherwise none skipped)')

    parser.add_argument('--n_threads', type=int, default=1, help='number of threads to run on')
    parser.add_argument('--timeout', type=int, default=600, help='task timeout in seconds')

    parser.add_argument('--max_meshes', type=int, default=-1)
    parser.add_argument('--use_ffield_cones', action='store_true', help='use cones from a *.ffield file (must have the same base name as the mesh file and be in dataset_dir)')
    parser.add_argument('--save_parameterized_meshes', action='store_true', help='save parameterized meshes as well as performance statistics. Meshes are placed in output_dir/meshes/')

    # Parse arguments
    args = parser.parse_args()

    mesh_dir = os.path.join(args.output_dir, 'meshes')
    ensure_dir_exists(args.output_dir)
    if args.save_parameterized_meshes:
        ensure_dir_exists(mesh_dir)

    # save the arguments
    with open(os.path.join(args.output_dir,"run_args.txt"), 'w') as f:
        argstr = " ".join(sys.argv)
        f.write(argstr)

    # Deal with lists
    if args.good_list:
        good_set = parse_file_list(args.good_list)
    if args.bad_list:
        bad_set = parse_file_list(args.bad_list)

    # Load the list of meshes
    meshes = []
    for f in os.listdir(args.dataset_dir):

        # respect lists
        f_name , f_ext = os.path.splitext(os.path.basename(f))
        if f_ext not in [".obj", ".stl", ".ply", ".off"]:
            continue

        if args.good_list and f_name not in good_set:
            continue
        if args.bad_list and f_name in bad_set:
            continue

        meshes.append(os.path.join(args.dataset_dir, f))
        if args.max_meshes > 0 and len(meshes) >= args.max_meshes:
            break

    print("Found {} input mesh files".format(len(meshes)))
    # random.shuffle(meshes)
    task_queue = ShellTaskQueue(nWorkers=args.n_threads, timeout=args.timeout)

    for m_path in meshes:
        m_basename = os.path.basename(m_path)
        m_base, m_ext = os.path.splitext(m_basename)

        if m_ext not in [".obj", ".stl", ".ply", ".off"]:
            continue

        output_path = os.path.join(args.output_dir, m_base + ".tsv")

        ffield_path = os.path.join(os.path.splitext(m_path)[0] + ".ffield")
        ffield_option = f"--ffield={ffield_path}" if args.use_ffield_cones else ""

        output_mesh_path = os.path.join(mesh_dir, m_base + ".obj")
        mesh_save_option = f"--outputMeshFilename={output_mesh_path}" if args.save_parameterized_meshes else ""

        cmd_list = [
            os.path.join(BIN_DIR,"parameterize"),
            m_path,
            ffield_option,
            mesh_save_option,
            f"--outputLogFilename={output_path}"
        ]

        # build the command
        cmd_str = " ".join(cmd_list)

        task_queue.add_task(cmd_str)

    task_queue.join()


if __name__ == "__main__":
    main()
