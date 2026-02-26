Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Parse-InvariantDouble {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Value
  )

  $parsed = 0.0
  $ok = [double]::TryParse(
    $Value,
    [System.Globalization.NumberStyles]::Any,
    [System.Globalization.CultureInfo]::InvariantCulture,
    [ref]$parsed
  )
  if (-not $ok) {
    throw "Cannot parse floating-point value '$Value' using InvariantCulture."
  }
  return $parsed
}

function New-RunResult {
  param(
    [Parameter(Mandatory = $true)][string]$RunId,
    [Parameter(Mandatory = $true)][int]$ExitCode,
    [Parameter(Mandatory = $true)][double]$ElapsedSeconds,
    [Parameter(Mandatory = $true)]$SummaryRow
  )

  if ($null -eq $SummaryRow) {
    return [PSCustomObject]@{
      run_id = $RunId
      status = "missing_summary"
      exit_code = $ExitCode
      elapsed_seconds = $ElapsedSeconds
      objective = [double]::PositiveInfinity
      total_travel_time = [double]::NaN
      budget = [double]::NaN
      budget_violation = [double]::NaN
      best_feasible_objective = [double]::NaN
      max_standard_iterations = 0
      max_optimality_iterations = 0
      route_search_threads = 0
      summary_timestamp = ""
    }
  }

  return [PSCustomObject]@{
    run_id = $RunId
    status = [string]$SummaryRow.Status
    exit_code = $ExitCode
    elapsed_seconds = $ElapsedSeconds
    objective = Parse-InvariantDouble -Value ([string]$SummaryRow.FinalObjective)
    total_travel_time = Parse-InvariantDouble -Value ([string]$SummaryRow.FinalTotalTravelTime)
    budget = Parse-InvariantDouble -Value ([string]$SummaryRow.FinalBudget)
    budget_violation = Parse-InvariantDouble -Value ([string]$SummaryRow.FinalBudgetViolation)
    best_feasible_objective = Parse-InvariantDouble -Value ([string]$SummaryRow.BestFeasibleObjective)
    max_standard_iterations = [int]$SummaryRow.MaxStandardIterations
    max_optimality_iterations = [int]$SummaryRow.MaxOptimalityIterations
    route_search_threads = [int]$SummaryRow.RouteSearchThreads
    summary_timestamp = [string]$SummaryRow.Timestamp
  }
}

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$exePath = Join-Path $projectRoot "build/mingw-vcpkg-release/main.exe"
$configPath = Join-Path $projectRoot "configs/cnd.siouxfalls.ini"
$summaryPath = Join-Path $projectRoot "performance_results/SiouxFalls/BilevelCND_run_summary.csv"
$optimalityPath = Join-Path $projectRoot "data/TransportationNetworks/SiouxFalls/SiouxFalls_optimality.csv"
$runOutputDir = Join-Path $projectRoot "performance_results/SiouxFalls/long_runs"
$etalonJsonPath = Join-Path $projectRoot "performance_results/SiouxFalls/etalon_solution_siouxfalls.json"
$etalonCsvPath = Join-Path $projectRoot "performance_results/SiouxFalls/etalon_runs_summary.csv"

if (-not (Test-Path $exePath)) {
  throw "Executable not found: $exePath"
}
if (-not (Test-Path $configPath)) {
  throw "Config not found: $configPath"
}

New-Item -ItemType Directory -Path $runOutputDir -Force | Out-Null

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$experiments = @(
  @{
    name = "A"
    run_id = "sf_ref_A_$stamp"
    args = @(
      "--max-standard-iters", "900",
      "--max-optcond-iters", "70",
      "--tolerance", "5e-5",
      "--link-threshold", "7e-4",
      "--budget-threshold", "8e-2",
      "--route-threads", "1"
    )
  },
  @{
    name = "B"
    run_id = "sf_ref_B_$stamp"
    args = @(
      "--max-standard-iters", "1100",
      "--max-optcond-iters", "60",
      "--tolerance", "1e-4",
      "--link-threshold", "5e-4",
      "--budget-threshold", "1e-1",
      "--route-threads", "1"
    )
  },
  @{
    name = "C"
    run_id = "sf_ref_C_$stamp"
    args = @(
      "--max-standard-iters", "800",
      "--max-optcond-iters", "90",
      "--tolerance", "3e-5",
      "--link-threshold", "1e-3",
      "--budget-threshold", "6e-2",
      "--route-threads", "1"
    )
  }
)

$results = @()

foreach ($exp in $experiments) {
  $runId = [string]$exp.run_id
  $stdoutPath = Join-Path $runOutputDir "$runId.stdout.log"
  $stderrPath = Join-Path $runOutputDir "$runId.stderr.log"
  Write-Host ""
  Write-Host ("[{0}] Starting run {1}" -f (Get-Date -Format "u"), $runId)
  Write-Host ("  Args: {0}" -f (($exp.args -join " ")))

  $sw = [System.Diagnostics.Stopwatch]::StartNew()
  & $exePath --config $configPath --metrics_run_id $runId --print-config false --loader-verbose false @($exp.args) 1> $stdoutPath 2> $stderrPath
  $exitCode = $LASTEXITCODE
  $sw.Stop()
  $elapsed = [Math]::Round($sw.Elapsed.TotalSeconds, 3)
  Write-Host ("[{0}] Finished run {1} exit={2} elapsed={3}s" -f (Get-Date -Format "u"), $runId, $exitCode, $elapsed)

  if (Test-Path $optimalityPath) {
    Copy-Item $optimalityPath (Join-Path $runOutputDir "$runId.optimality.csv") -Force
  }

  $summaryRow = $null
  if (Test-Path $summaryPath) {
    $summaryRow = Import-Csv $summaryPath | Where-Object { $_.RunId -eq $runId } | Select-Object -Last 1
  }
  $results += New-RunResult -RunId $runId -ExitCode $exitCode -ElapsedSeconds $elapsed -SummaryRow $summaryRow
}

$successResults = $results |
  Where-Object { $_.status -eq "success" -and $_.exit_code -eq 0 -and [double]::IsFinite([double]$_.objective) }

if ($successResults.Count -eq 0) {
  $results | Export-Csv -Path $etalonCsvPath -NoTypeInformation -Encoding UTF8
  throw "No successful runs with finite objective. See logs in $runOutputDir"
}

$best = $successResults | Sort-Object objective | Select-Object -First 1
$bestOptimalityPath = Join-Path $runOutputDir ("{0}.optimality.csv" -f $best.run_id)

if (-not (Test-Path $bestOptimalityPath)) {
  throw "Best run optimality file not found: $bestOptimalityPath"
}

$optimalityRows = Import-Csv $bestOptimalityPath
$capacities = @()
for ($i = 0; $i -lt $optimalityRows.Count; ++$i) {
  $row = $optimalityRows[$i]
  $capacities += [PSCustomObject]@{
    link_index = $i
    flow = Parse-InvariantDouble -Value ([string]$row.flow)
    lower_bound = Parse-InvariantDouble -Value ([string]$row.lower_bound)
    upper_bound = Parse-InvariantDouble -Value ([string]$row.upper_bound)
    capacity = Parse-InvariantDouble -Value ([string]$row.capacity)
    condition_result = Parse-InvariantDouble -Value ([string]$row.condition_result)
  }
}

$etalon = [PSCustomObject]@{
  generated_at = (Get-Date).ToString("o")
  dataset = "SiouxFalls"
  etalon_run_id = $best.run_id
  etalon_objective = [double]$best.objective
  etalon_total_travel_time = [double]$best.total_travel_time
  etalon_budget = [double]$best.budget
  etalon_best_feasible_objective = [double]$best.best_feasible_objective
  etalon_elapsed_seconds = [double]$best.elapsed_seconds
  route_search_threads = [int]$best.route_search_threads
  max_standard_iterations = [int]$best.max_standard_iterations
  max_optimality_iterations = [int]$best.max_optimality_iterations
  run_output_dir = $runOutputDir
  runs = $results
  capacities = $capacities
}

$etalon | ConvertTo-Json -Depth 8 | Out-File -FilePath $etalonJsonPath -Encoding UTF8
$results | Export-Csv -Path $etalonCsvPath -NoTypeInformation -Encoding UTF8

Write-Host ""
Write-Host "Completed long-run benchmark."
Write-Host ("Etalon run: {0}" -f $best.run_id)
Write-Host ("Etalon objective: {0}" -f $best.objective)
Write-Host ("Saved: {0}" -f $etalonJsonPath)
Write-Host ("Saved: {0}" -f $etalonCsvPath)
