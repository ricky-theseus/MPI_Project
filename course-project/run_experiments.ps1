param([int]$runs=10, [int]$gen=200000, [int]$np=4)

$start = Get-Date
Write-Host "TSP Experiments: $runs runs x $gen gens x $np procs"
Write-Host "Start: $start"

$algs = @("dEA", "HdEA", "HdEA_MP")
$dirs = @(
    "01_分布式进化算法",
    "02_分层分布式进化算法",
    "03_移动种群分层分布式进化算法"
)

for ($i = 0; $i -lt $algs.Length; $i++) {
    $a = $algs[$i]
    $d = $dirs[$i]
    $exe = "$d\TSP_$a.exe"
    $res = "results_$a.txt"
    $resPath = Join-Path (Get-Location) $res
    if (Test-Path $resPath) { Remove-Item $resPath }

    Write-Host "--- Running $a ---"
    $t0 = Get-Date
    $p = Start-Process mpiexec -ArgumentList "-n $np `"$exe`" $runs $gen" -NoNewWindow -Wait -PassThru
    $t1 = Get-Date
    $dur = ($t1 - $t0).TotalMinutes
    Write-Host "--- $a done in $dur min (exit $($p.ExitCode)) ---"
}

$end = Get-Date
Write-Host "Total: $(($end - $start).TotalMinutes) min"
