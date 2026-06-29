#!/usr/bin/env python3
import sys
import os
import time


def get_progress(log_file):
    if not os.path.exists(log_file):
        return 0, False
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            done_count = content.count("completata:")
            is_finished = ("Simulation completed" in content or "Optimization completed" in content
                           or "DoF scan completed" in content or "PSF+DoF scan completed" in content)
            return done_count, is_finished
    except Exception:
        return 0, False


def main_tty(tmpdir, n_jobs, totals):
    """Dashboard rich per uso interattivo (stdout è un terminale)."""
    from rich.progress import (Progress, BarColumn, TextColumn,
                                TimeElapsedColumn, TimeRemainingColumn,
                                SpinnerColumn, MofNCompleteColumn)
    from rich.console import Console

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
        transient=False,
    ) as progress:
        job_tasks = []
        for i in range(n_jobs):
            task_id = progress.add_task(f"Chunk {i:02d}", total=totals[i] if i < len(totals) else 1)
            job_tasks.append(task_id)

        overall_total = sum(totals)
        overall_task = progress.add_task("[bold green]Overall Progress",
                                         total=overall_total if overall_total > 0 else 1)

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
            time.sleep(1)

    console.print("\n[bold green]Tutti i job completati![/bold green]")


def main_pipe(tmpdir, n_jobs, totals):
    """Progress plain-text per uso da lens_runner.py (stdout non è tty)."""
    total = sum(totals)
    t_start = time.time()
    last_pct = -1
    last_emit_elapsed = 0.0

    while True:
        all_done = True
        done = 0
        for i in range(n_jobs):
            d, finished = get_progress(os.path.join(tmpdir, f"chunk_{i}.log"))
            done += d
            if not finished:
                all_done = False

        elapsed = time.time() - t_start
        pct = int(done * 100 / total) if total > 0 else 0
        eta = elapsed * (total - done) / done if done > 0 else 0.0

        should_emit = (
            pct >= last_pct + 10
            or (elapsed - last_emit_elapsed >= 300 and done > 0 and pct > last_pct)
        )
        if should_emit:
            print(f"[PROGRESS] {pct}% ({done}/{total} run) — elapsed {elapsed:.0f}s ETA {eta:.0f}s",
                  flush=True)
            last_pct = (pct // 10) * 10
            last_emit_elapsed = elapsed

        if all_done:
            print(f"[PROGRESS] 100% ({done}/{total} run) — completato in {elapsed:.0f}s",
                  flush=True)
            break

        time.sleep(30)


def main():
    if len(sys.argv) < 4:
        sys.exit(1)

    tmpdir = sys.argv[1]
    n_jobs = int(sys.argv[2])
    totals = [int(x) for x in sys.argv[3:]]

    if sys.stdout.isatty():
        main_tty(tmpdir, n_jobs, totals)
    else:
        main_pipe(tmpdir, n_jobs, totals)


if __name__ == "__main__":
    main()
