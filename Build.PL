use 5.006;
use strict;
use warnings;
use Module::Build;

# loading bin scripts #
opendir IN, "./bin/" or die $!;
my @binfiles = readdir IN;
my @binscripts = grep(/.pl$|.py$|.r$/i, @binfiles);
@binscripts = map{$_ = "./bin/" . $_} @binscripts;

my $builder = Module::Build->new(
    module_name         => 'deltaGC',
    license             => 'perl',
    dist_author         => q{Nick Youngblut <ndy2@cornell.edu>},
    dist_version_from   => 'lib/deltaGC.pm',
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        'perl' => 5.006,
    },
    add_to_cleanup      => [ 'deltaGC-*' ],
    create_makefile_pl => 'traditional',
    script_files => \@binscripts			   
);

$builder->create_build_script();
