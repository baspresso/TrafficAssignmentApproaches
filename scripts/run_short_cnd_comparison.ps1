param(
  [string]$EtalonRunId = "1772049476220"
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

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$exePath = Join-Path $projectRoot "build/mingw-vcpkg-release/main.exe"
$configPath = Join-Path $projectRoot "configs/cnd.siouxfalls.ini"
$summaryPath = Join-Path $projectRoot "performance_results/SiouxFalls/BilevelCND_run_summary.csv"
$optSourcePath = Join-Path $projectRoot "data/TransportationNetworks/SiouxFalls/SiouxFalls_optimality.csv"
$outDir = Join-Path $projectRoot "performance_results/SiouxFalls/short_runs"

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

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$experiments = @(
  @{
    scenario = "optcond_only"
    run_id = "sf_short_optcond_$stamp"
    args = @(
      "--nlopt-algorithm", "LN_COBYLA",
      "--max-standard-iters", "0",
      "--max-optcond-iters", "25",
      "--route-threads", "1",
      "--tolerance", "1e-4"
    )
  },
  @{
    scenario = "standard_only_LN_COBYLA"
    run_id = "sf_short_std_cobyla_$stamp"
    args = @(
      "--nlopt-algorithm", "LN_COBYLA",
      "--max-standard-iters", "220",
      "--max-optcond-iters", "0",
      "--route-threads", "1",
      "--tolerance", "1e-4"
    )
  },
  @{
    scenario = "standard_only_GN_ISRES"
    run_id = "sf_short_std_isres_$stamp"
    args = @(
      "--nlopt-algorithm", "GN_ISRES",
      "--max-standard-iters", "220",
      "--max-optcond-iters", "0",
      "--route-threads", "1",
      "--tolerance", "1e-4"
    )
  },
  @{
    scenario = "standard_only_AUGLAG"
    run_id = "sf_short_std_auglag_$stamp"
    args = @(
      "--nlopt-algorithm", "AUGLAG",
      "--max-standard-iters", "220",
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

  $sw = [System.Diagnostics.Stopwatch]::StartNew()
  & $exePath --config $configPath --metrics_run_id $runId --print-config false --loader-verbose false @($exp.args) 1> $stdoutPath 2> $stderrPath
  $exitCode = $LASTEXITCODE
  $sw.Stop()
  $elapsed = [Math]::Round($sw.Elapsed.TotalSeconds, 3)
  Write-Host ("Finished {0}; exit={1}; elapsed={2}s" -f $runId, $exitCode, $elapsed)

  if (Test-Path $optSourcePath) {
    Copy-Item $optSourcePath (Join-Path $outDir "$runId.optimality.csv") -Force
  }

  $row = (Import-Csv $summaryPath | Where-Object { $_.RunId -eq $runId } | Select-Object -Last 1)
  if ($null -eq $row) {
    $runResults += [PSCustomObject]@{
      scenario = $exp.scenario
      run_id = $runId
      status = "missing_summary"
      exit_code = $exitCode
      elapsed_seconds = $elapsed
      objective = [double]::NaN
      objective_delta_to_etalon = [double]::NaN
      total_travel_time = [double]::NaN
      budget = [double]::NaN
      algorithm = ""
      max_standard_iterations = 0
      max_optimality_iterations = 0
      route_search_threads = 0
    }
    continue
  }

  $objective = Parse-InvariantDouble ([string]$row.FinalObjective)
  $runResults += [PSCustomObject]@{
    scenario = $exp.scenario
    run_id = $runId
    status = [string]$row.Status
    exit_code = $exitCode
    elapsed_seconds = $elapsed
    objective = $objective
    objective_delta_to_etalon = $objective - $etalonObjective
    total_travel_time = Parse-InvariantDouble ([string]$row.FinalTotalTravelTime)
    budget = Parse-InvariantDouble ([string]$row.FinalBudget)
    algorithm = [string]$row.Algorithm
    max_standard_iterations = [int]$row.MaxStandardIterations
    max_optimality_iterations = [int]$row.MaxOptimalityIterations
    route_search_threads = [int]$row.RouteSearchThreads
  }
}

$resultCsvPath = Join-Path $outDir ("short_comparison_runs_{0}.csv" -f $stamp)
$runResults | Export-Csv -Path $resultCsvPath -NoTypeInformation -Encoding UTF8

Write-Host ""
Write-Host "Completed short comparison runs."
Write-Host ("Etalon objective: {0}" -f $etalonObjective)
Write-Host ("Result CSV: {0}" -f $resultCsvPath)
$runResults | Format-Table -AutoSize
