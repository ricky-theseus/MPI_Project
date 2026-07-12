param(
    [int]$Runs = 10,
    [long]$Gens = 200000
)

$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ResultsDir = "$ProjectDir\results"

@("serial", "dEA", "HdEA", "HdEA_MP") | ForEach-Object {
    $dir = "$ResultsDir\$_"
    if (-not (Test-Path $dir)) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }
    $convDir = "$dir\convergence"
    if (-not (Test-Path $convDir)) { New-Item -ItemType Directory -Path $convDir -Force | Out-Null }
}

Write-Host "=== Running experiments: $Runs runs, $Gens generations ===" -ForegroundColor Cyan

$MpiOpts = @("-n", "4", "-wdir", $ProjectDir)

# Step 1: Serial baseline
Write-Host "`n[1/4] Serial baseline (TSP0)..." -ForegroundColor Yellow
if (Test-Path "$ProjectDir\TSP0.exe") {
    for ($r = 1; $r -le $Runs; $r++) {
        Write-Host "  Run $r/$Runs..."
        Push-Location $ProjectDir
        & "$ProjectDir\TSP0.exe"
        Pop-Location
    }
} else {
    Write-Host "  WARNING: TSP0.exe not found, skipping serial baseline" -ForegroundColor Red
}

# Step 2: dEA
Write-Host "`n[2/4] dEA (4 processes)..." -ForegroundColor Yellow
& mpiexec @MpiOpts "$ProjectDir\TSP_dEA.exe" $Runs $Gens
$deaIr = "$ResultsDir\dEA"
if (Test-Path "$ProjectDir\results_dEA.txt") { Move-Item -Force "$ProjectDir\results_dEA.txt" "$deaIr\results.txt" }
if (Test-Path "$ProjectDir\convergence_dEA.txt") { Move-Item -Force "$ProjectDir\convergence_dEA.txt" "$deaIr\convergence.txt" }

# Step 3: HdEA
Write-Host "`n[3/4] HdEA (4 processes, 4 islands)..." -ForegroundColor Yellow
& mpiexec @MpiOpts "$ProjectDir\TSP_HdEA.exe" $Runs $Gens
$hdeaDir = "$ResultsDir\HdEA"
if (Test-Path "$ProjectDir\results_HdEA.txt") { Move-Item -Force "$ProjectDir\results_HdEA.txt" "$hdeaDir\results.txt" }
if (Test-Path "$ProjectDir\convergence_HdEA.txt") { Move-Item -Force "$ProjectDir\convergence_HdEA.txt" "$hdeaDir\convergence.txt" }

# Step 4: HdEA-MP
Write-Host "`n[4/4] HdEA-MP (4 processes, moving colony)..." -ForegroundColor Yellow
& mpiexec @MpiOpts "$ProjectDir\TSP_HdEA_MP.exe" $Runs $Gens
$hdeampDir = "$ResultsDir\HdEA_MP"
if (Test-Path "$ProjectDir\results_HdEA_MP.txt") { Move-Item -Force "$ProjectDir\results_HdEA_MP.txt" "$hdeampDir\results.txt" }
if (Test-Path "$ProjectDir\convergence_HdEA_MP.txt") { Move-Item -Force "$ProjectDir\convergence_HdEA_MP.txt" "$hdeampDir\convergence.txt" }

Write-Host "`n=== All experiments complete ===" -ForegroundColor Cyan
Write-Host "Results saved to $ResultsDir"
