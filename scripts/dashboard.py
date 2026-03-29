#!/usr/bin/env python3
import sys
import os
import time
import re
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn, SpinnerColumn, MofNCompleteColumn
from rich.console import Console
from rich.table import Table
from rich.live import Live
from rich.panel import Panel
from rich.layout import Layout

def get_progress(log_file):
    if not os.path.exists(log_file):
        return 0, False
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            done_count = content.count("Run done:")
            is_finished = "Simulation completed" in content or "Optimization completed" in content
            return done_count, is_finished
    except:
        return 0, False

def main():
    if len(sys.argv) < 4:
        sys.exit(1)

    tmpdir = sys.argv[1]
    n_jobs = int(sys.argv[2])
    # The rest are totals
    totals = [int(x) for x in sys.argv[3:]]

    console = Console()
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
        transient=False
    ) as progress:
        
        job_tasks = []
        for i in range(n_jobs):
            task_id = progress.add_task(f"Chunk {i:02d}", total=totals[i] if i < len(totals) else 1)
            job_tasks.append(task_id)
        
        overall_total = sum(totals)
        overall_task = progress.add_task("[bold green]Overall Progress", total=overall_total if overall_total > 0 else 1)

        start_time = time.time()
        while True:
            all_done = True
            total_done = 0
            
            for i in range(n_jobs):
                log_file = os.path.join(tmpdir, f"chunk_{i}.log")
                done, finished = get_progress(log_file)
                
                progress.update(job_tasks[i], completed=done)
                total_done += done
                
                if not finished:
                    all_done = False
            
            progress.update(overall_task, completed=total_done)
            
            if all_done:
                break
                
            # Check if all processes are still alive (optional, but good for error detection)
            # For now, just sleep
            time.sleep(1)

    console.print("\n[bold green]✅ All jobs completed successfully![/bold green]")

if __name__ == "__main__":
    main()
