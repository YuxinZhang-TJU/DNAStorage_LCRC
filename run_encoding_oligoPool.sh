#!/bin/bash

export PATH=bin:$PATH

out_dir=output

if [ ! -d $out_dir  ]; then
	mkdir -p $out_dir
fi

Encode_OligoPool $out_dir/oligo_pool.fa
