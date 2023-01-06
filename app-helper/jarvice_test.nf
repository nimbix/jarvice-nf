process stage1 {
  debug true
  executor 'jarvice'
  container 'us-docker.pkg.dev/jarvice/images/ubuntu-desktop:bionic'
  machineType 'micro'
  cpus 1

  input:
      val x

  output:
    path 'stage1_*.txt'

  """
  sleep 5
  i=\$RANDOM
  echo "Processed $x \$i" > stage1_\$i.txt
  """
}


workflow {
  Channel.from(1..2) | map { x -> "input_$x"} | stage1 | view { it }
}
