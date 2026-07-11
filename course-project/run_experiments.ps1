param(
    [int]$Runs = 10,
    [long]$Gens = 200000
)

$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ResultsDir = "$ProjectDir\results"

# Create results directories
@("serial", "dEA", "HdEA", "HdEA_MP") | ForEach-Object {
    $dir = "$ResultsDir\$_"
    if (-not (Test-Path $dir)) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }
    $convDir = "$dir\convergence"
    if (-not (Test-Path $convDir)) { New-Item -ItemType Directory -Path $convDir -Force | Out-Null }
}

Write-Host "=== Running experiments: $Runs runs, $Gens generations ===" -ForegroundColor Cyan

# Step 1: Serial baseline (TSP0.C compiled as TSP0.exe)
Write-Host "`n[1/4] Serial baseline (TSP0)..." -ForegroundColor Yellow
$serialExe = "$ProjectDir\TSP0.exe"
if (Test-Path $serialExe) {
    for ($r = 1; $r -le $Runs; $r++) {
        Write-Host "  Run $r/$Runs..."
        $outFile = "$ResultsDir\serial\convergence\run_$($r.ToString('00')).txt"
        & $serialExe $r $Gens
        # TSP0 appends to e:\tsp0.txt, so we need its output
    }
} else {
    Write-Host "  WARNING: TSP0.exe not found, skipping serial baseline" -ForegroundColor Red
}

# Step 2: dEA
Write-Host "`n[2/4] dEA (4 processes)..." -ForegroundColor Yellow
$deaIr = "$ResultsDir\dEA"
& mpiexec -n 4 "$ProjectDir\01_分布式进化算法\TSP_dEA.exe" $Runs $Gens
if (Test-Path "$ProjectDir\01_分布式进化算法\results_dEA.txt") {
    Move-Item -Force "$ProjectDir\01_分布式进化算法\results_dEA.txt" "$deaIr\results.txt"
}
if (Test-Path "$ProjectDir\01_分布式进化算法\convergence_dEA.txt") {
    Move-Item -Force "$ProjectDir\01_分布式进化算法\convergence_dEA.txt" "$deaIr\convergence.txt"
}

# Step 3: HdEA
Write-Host "`n[3/4] HdEA (4 processes, 4 islands)..." -ForegroundColor Yellow
$hdeaDir = "$ResultsDir\HdEA"
& mpiexec -n 4 "$ProjectDir\02_分层分布式进化算法\TSP_HdEA.exe" $Runs $Gens
if (Test-Path "$ProjectDir\02_分层分布式进化算法\results_HdEA.txt") {
    Move-Item -Force "$ProjectDir\02_分层分布式进化算法\results_HdEA.txt" "$hdeaDir\results.txt"
}
if (Test-Path "$ProjectDir\02_分层分布式进化算法\convergence_HdEA.txt") {
    Move-Item -Force "$ProjectDir\02_分层分布式进化算法\convergence_HdEA.txt" "$hdeaDir\convergence.txt"
}

# Step 4: HdEA-MP
Write-Host "`n[4/4] HdEA-MP (4 processes, moving colony)..." -ForegroundColor Yellow
$hdeampDir = "$ResultsDir\HdEA_MP"
& mpiexec -n 4 "$ProjectDir\03_HdEA_MP\TSP_HdEA_MP.exe" $Runs $Gens
if (Test-Path "$ProjectDir\03_HdEA_MP\results_HdEA_MP.txt") {
    Move-Item -Force "$ProjectDir\03_HdEA_MP\results_HdEA_MP.txt" "$hdeampDir\results.txt"
}
if (Test-Path "$ProjectDir\03_HdEA_MP\convergence_HdEA_MP.txt") {
    Move-Item -Force "$ProjectDir\03_HdEA_MP\convergence_HdEA_MP.txt" "$hdeampDir\convergence.txt"
}

Write-Host "`n=== All experiments complete ===" -ForegroundColor Cyan
Write-Host "Results saved to $ResultsDir"
