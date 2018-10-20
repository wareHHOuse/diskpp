#!/usr/bin/perl
use warnings;
use strict;

my @bi_nums = (2, 50);
my @alphas = (0.1, 1, 10, 100);
my @meshes = (32);

sub prepare_sim {
    my $bi_num  = shift;
    my $alpha   = shift;
    my $mesh_h  = shift;
    
    my $bingham_config = <<"    EOF";
config.degree_cell = 0
config.degree_face = 0
    
config.input_mesh = "../../../../diskpp/meshes/2D_quads/diskpp/testmesh-$mesh_h-$mesh_h.quad"
    
bi.hname = "$mesh_h";
bi.alpha = $alpha;  --ALG augmentation parameter
bi.Lref  = 1;       --Reference dimension
bi.Vref  = 1;       --Reference velocity
bi.mu = 1;          --Viscosity
bi.Bn = $bi_num;    --Bingham number
bi.f  = 0;          --force
bi.problem = "DRIVEN" --Choose only "circular","annulus", or "square"
    EOF
    
    my $dirname  = "bingham_" . $bi_num . "_" . $alpha . "_". $mesh_h;
    my $filename = "$dirname/config.lua";

    mkdir $dirname;
    
    open(CONFIG, ">$filename");
    print CONFIG $bingham_config;
    link "../bingham/bingham_vector", "$dirname/bingham_vector";
}

foreach my $bi_num (@bi_nums) {
    foreach my $alpha (@alphas) {
        foreach my $mesh (@meshes) {
            prepare_sim($bi_num, $alpha, $mesh);
        }
    }
}


