param(
  [string]$EtalonRunId = "1772049476220",
  [int]$PerRunTimeoutSeconds = 300
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Parse-InvariantDouble {
  param([Parameter(Mandatory = $true)][string]$Value)
  $result = 0.0
  if (-not [double]::TryParse(
      $Value,
      [System.Globalization.NumberStyles]::Any,
      [System.Globalization.CultureInfo]::InvariantCulture,
      [ref]$result
    )) {
    throw "Cannot parse floating-point value '$Value'."
  }
  return $result
}

function TryParse-InvariantDouble {
  param([string]$Value)
  if ([string]::IsNullOrWhiteSpace($Value)) {
    return [double]::NaN
  }
  $result = 0.0
  if ([double]::TryParse(
      $Value,
      [System.Globalization.NumberStyles]::Any,
      [System.Globalization.CultureInfo]::InvariantCulture,
      [ref]$result
    )) {
    return $result
  }
  return [double]::NaN
}

function Get-RunSummaryRow {
  param(
    [Parameter(Mandatory = $true)][string]$SummaryPath,
    [Parameter(Mandatory = $true)][string]$RunId
  )
  if (-not (Test-Path $SummaryPath)) {
    return $null
  }
  return (Import-Csv $SummaryPath | Where-Object { $_.RunId -eq $RunId } | Select-Object -Last 1)
}

function Get-TraceSnapshot {
  param(
    [Parameter(Mandatory = $true)][string]$TraceRoot,
    [Parameter(Mandatory = $true)][string]$RunId
  )

  $traceFile = Get-ChildItem -Path $TraceRoot -Filter "*_${RunId}_quality_time.csv" -File -ErrorAction SilentlyContinue |
    Select-Object -First 1
  if ($null -eq $traceFile) {
    return [PSCustomObject]@{
      has_trace = $false
      trace_file = ""
      best_feasible_objective = [double]::NaN
      last_objective = [double]::NaN
      last_total_travel_time = [double]::NaN
      last_budget = [double]::NaN
      last_elapsed_seconds = [double]::NaN
    }
  }

  $rows = @(Import-Csv $traceFile.FullName)
  if ($rows.Count -eq 0) {
    return [PSCustomObject]@{
      has_trace = $true
      trace_file = $traceFile.FullName
      best_feasible_objective = [double]::NaN
      last_objective = [double]::NaN
      last_total_travel_time = [double]::NaN
      last_budget = [double]::NaN
      last_elapsed_seconds = [double]::NaN
    }
  }

  $bestFeasible = [double]::NaN
  foreach ($r in $rows) {
    $objective = TryParse-InvariantDouble ([string]$r.Objective)
    $budgetViolation = TryParse-InvariantDouble ([string]$r.BudgetViolation)
    if ([double]::IsFinite($objective) -and [double]::IsFinite($budgetViolation) -and $budgetViolation -le 1e-10) {
      if (-not [double]::IsFinite($bestFeasible) -or $objective -lt $bestFeasible) {
        $bestFeasible = $objective
      }
    }
  }

  $last = $rows | Select-Object -Last 1
  return [PSCustomObject]@{
    has_trace = $true
    trace_file = $traceFile.FullName
    best_feasible_objective = $bestFeasible
    last_objective = TryParse-InvariantDouble ([string]$last.Objective)
    last_total_travel_time = TryParse-InvariantDouble ([string]$last.TotalTravelTime)
    last_budget = TryParse-InvariantDouble ([string]$last.Budget)
    last_elapsed_seconds = TryParse-InvariantDouble ([string]$last.'ElapsedTime(s)')
  }
}

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$exePath = Join-Path $projectRoot "build/mingw-vcpkg-release/main.exe"
$configPath = Join-Path $projectRoot "configs/cnd.siouxfalls.ini"
$summaryPath = Join-Path $projectRoot "performance_results/SiouxFalls/BilevelCND_run_summary.csv"
$traceRoot = Join-Path $projectRoot "performance_results/SiouxFalls"
$optSourcePath = Join-Path $projectRoot "data/TransportationNetworks/SiouxFalls/SiouxFalls_optimality.csv"
$outDir = Join-Path $projectRoot "performance_results/SiouxFalls/five_min_runs"

if (-not (Test-Path $exePath)) { throw "Executable not found: $exePath" }
if (-not (Test-Path $configPath)) { throw "Config not found: $configPath" }
if (-not (Test-Path $summaryPath)) { throw "Summary CSV not found: $summaryPath" }

New-Item -ItemType Directory -Path $outDir -Force | Out-Null

$summaryRows = Import-Csv $summaryPath
$etalonRow = $summaryRows | Where-Object { $_.RunId -eq $EtalonRunId } | Select-Object -Last 1
if ($null -eq $etalonRow) {
  throw "Etalon run_id '$EtalonRunId' not found in $summaryPath"
}

$etalonObjective = Parse-InvariantDouble ([string]$etalonRow.FinalObjective)
Write-Host ("Etalon run_id={0}, objective={1}" -f $EtalonRunId, $etalonObjective)
Write-Host ("Per-run timeout: {0}s" -f $PerRunTimeoutSeconds)

# Ensure no stale process from previous interrupted run.
Get-Process main -ErrorAction SilentlyContinue | Stop-Process -Force

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$experiments = @(
  @{
    scenario = "optcond_only"
    algorithm_name = "LN_COBYLA"
    run_id = "sf_5min_optcond_$stamp"
    args = @(
      "--nlopt-algorithm", "LN_COBYLA",
      "--max-standard-iters", "0",
      "--max-optcond-iters", "5000",
      "--route-threads", "1",
      "--tolerance", "1e-4"
    )
  },
  @{
    scenario = "standard_only_LN_COBYLA"
    algorithm_name = "LN_COBYLA"
    run_id = "sf_5min_std_cobyla_$stamp"
    args = @(
      "--nlopt-algorithm", "LN_COBYLA",
      "--max-standard-iters", "5000",
      "--max-optcond-iters", "0",
      "--route-threads", "1",
      "--tolerance", "1e-4"
    )
  },
  @{
    scenario = "standard_only_LN_NELDERMEAD"
    algorithm_name = "LN_NELDERMEAD"
    run_id = "sf_5min_std_neldermead_$stamp"
    args = @(
      "--nlopt-algorithm", "LN_NELDERMEAD",
      "--max-standard-iters", "5000",
      "--max-optcond-iters", "0",
      "--route-threads", "1",
      "--tolerance", "1e-4"
    )
  },
  @{
    scenario = "standard_only_LN_SBPLX"
    algorithm_name = "LN_SBPLX"
    run_id = "sf_5min_std_sbplx_$stamp"
    args = @(
      "--nlopt-algorithm", "LN_SBPLX",
      "--max-standard-iters", "5000",
      "--max-optcond-iters", "0",
      "--route-threads", "1",
      "--tolerance", "1e-4"
    )
  }
)

$runResults = @()

foreach ($exp in $experiments) {
  $runId = [string]$exp.run_id
  $stdoutPath = Join-Path $outDir "$runId.stdout.log"
  $stderrPath = Join-Path $outDir "$runId.stderr.log"
  Write-Host ""
  Write-Host ("Starting {0} ({1})" -f $runId, $exp.scenario)

  $allArgs = @(
    "--config", $configPath,
    "--metrics_run_id", $runId,
    "--metrics_flush_every_n_points", "1",
    "--print-config", "false",
    "--loader-verbose", "false"
  ) + @($exp.args)
  $sw = [System.Diagnostics.Stopwatch]::StartNew()
  $proc = Start-Process -FilePath $exePath `
                        -ArgumentList $allArgs `
                        -RedirectStandardOutput $stdoutPath `
                        -RedirectStandardError $stderrPath `
                        -NoNewWindow `
                        -PassThru

  $completed = $proc.WaitForExit($PerRunTimeoutSeconds * 1000)
  $timedOut = -not $completed
  if ($timedOut) {
    & taskkill /PID $proc.Id /T /F | Out-Null
    Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
    Get-Process main -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue

    $deadline = (Get-Date).AddSeconds(10)
    while ((Get-Process -Id $proc.Id -ErrorAction SilentlyContinue) -and (Get-Date) -lt $deadline) {
      Start-Sleep -Milliseconds 200
    }
  }
  $sw.Stop()

  $exitCode = if ($timedOut) { 124 } else { $proc.ExitCode }
  $elapsed = [Math]::Round($sw.Elapsed.TotalSeconds, 3)
  Write-Host ("Finished {0}; exit={1}; elapsed={2}s; timeout={3}" -f $runId, $exitCode, $elapsed, $timedOut)

  if (Test-Path $optSourcePath) {
    Copy-Item $optSourcePath (Join-Path $outDir "$runId.optimality.csv") -Force
  }

  $summaryRow = Get-RunSummaryRow -SummaryPath $summaryPath -RunId $runId
  $traceInfo = Get-TraceSnapshot -TraceRoot $traceRoot -RunId $runId

  $status = ""
  $objectiveForComparison = [double]::NaN
  $objective = [double]::NaN
  $travelTime = [double]::NaN
  $budget = [double]::NaN
  $algorithm = [string]$exp.algorithm_name
  $maxStandard = 0
  $maxOpt = 0
  $routeThreads = 0

  if ($null -ne $summaryRow) {
    $status = [string]$summaryRow.Status
    if ($timedOut) {
      $status = "timeout_with_partial_summary"
    }
    $objective = TryParse-InvariantDouble ([string]$summaryRow.FinalObjective)
    $travelTime = TryParse-InvariantDouble ([string]$summaryRow.FinalTotalTravelTime)
    $budget = TryParse-InvariantDouble ([string]$summaryRow.FinalBudget)
    $algorithm = [string]$summaryRow.Algorithm
    $maxStandard = [int]$summaryRow.MaxStandardIterations
    $maxOpt = [int]$summaryRow.MaxOptimalityIterations
    $routeThreads = [int]$summaryRow.RouteSearchThreads
  } else {
    $status = if ($timedOut) { "timeout_missing_summary" } else { "missing_summary" }
    $objective = $traceInfo.last_objective
    $travelTime = $traceInfo.last_total_travel_time
    $budget = $traceInfo.last_budget
  }

  if ($traceInfo.has_trace -and [double]::IsFinite($traceInfo.best_feasible_objective)) {
    $objectiveForComparison = $traceInfo.best_feasible_objective
  } else {
    $objectiveForComparison = $objective
  }

  $objectiveDelta = [double]::NaN
  if ([double]::IsFinite($objectiveForComparison)) {
    $objectiveDelta = $objectiveForComparison - $etalonObjective
  }

  $runResults += [PSCustomObject]@{
    scenario = $exp.scenario
    run_id = $runId
    status = $status
    timed_out = $timedOut
    exit_code = $exitCode
    elapsed_seconds = $elapsed
    objective = $objective
    objective_for_comparison = $objectiveForComparison
    objective_delta_to_etalon = $objectiveDelta
    total_travel_time = $travelTime
    budget = $budget
    trace_best_feasible_objective = $traceInfo.best_feasible_objective
    trace_last_elapsed_seconds = $traceInfo.last_elapsed_seconds
    trace_file = $traceInfo.trace_file
    algorithm = $algorithm
    max_standard_iterations = $maxStandard
    max_optimality_iterations = $maxOpt
    route_search_threads = $routeThreads
  }
}

$resultCsvPath = Join-Path $outDir ("five_min_comparison_runs_{0}.csv" -f $stamp)
$runResults | Export-Csv -Path $resultCsvPath -NoTypeInformation -Encoding UTF8

Write-Host ""
Write-Host "Completed 5-minute comparison runs."
Write-Host ("Etalon objective: {0}" -f $etalonObjective)
Write-Host ("Result CSV: {0}" -f $resultCsvPath)
$runResults | Format-Table -AutoSize
