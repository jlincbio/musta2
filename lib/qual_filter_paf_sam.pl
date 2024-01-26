#!/usr/bin/env perl

## quality_filter_paf_sam.pl
## 4/26/2022

use strict;
use warnings;
use Sort::Key::Natural 'natsort';

my $prog_id = uc 'QUAL FILTERING';

sub get_index {
	my ($header, @ids) = @_;
	my @x = @{$header};
	my %k;
	foreach my $y (@ids) {
		my @z = grep {$x[$_] eq $y } 0..$#x;
		if ($#z > 0) {
			@z = sort {$a <=> $b} @z;
		}
		$k{$y} = (shift @z);
	}
	return \%k;
}
sub read_paf_quality {
	my ($f, $q) = @_;
	open my $F_INPUT, '<', $f or die "### $prog_id ERROR: can't open $f\n";
	my (%map_quality, %p);
	while (<$F_INPUT>) {
		chomp;
		my $p_cur = $_;
		my @s = split '\t', $p_cur;
		# seq_id -> 0, map_quality -> 11 
	    # c__________i____________
	    # query_start, query_end, start, end is 0-based
		push @{$map_quality{$s[0]}}, $s[11];
		push @{$p{$s[0]}}, $p_cur;
	}
	close $F_INPUT;
	
	my @ids;
	foreach my $id (keys %map_quality) {
		my $n_pass = 0;
		my @mapq = @{$map_quality{$id}};
		foreach my $q_cur (@mapq) {
			$n_pass++ if ($q_cur > $q);
		}
		if ($n_pass >= scalar @mapq) {
			push @ids, $id;
		}
	}
	my %p1;
	foreach my $id (@ids) {
		$p1{$id} = $p{$id};
	}
	return (\%map_quality, \@ids, \%p1);
}
sub qfilt {
	my ($paf, $sam, $out_paf, $out_sam, $thres, $s_name) = @_;
	my ($r_paf_quality, $r_pass_id, $r_paf) = read_paf_quality($paf, $thres);
	my %pass_id = map {$_ => 1} @{$r_pass_id};
	my @pass_id_sorted = natsort (keys %pass_id);
	
	my %pass_paf = %{$r_paf};
	open my $F_PAF_OUTPUT, '>', $out_paf or die "### $prog_id ERROR: can't write filtered PAF for sample $s_name!\n";
	foreach my $id (@pass_id_sorted) {
		#my $q_id = $r_paf_quality->{$id};
		#foreach my $y (@{$q_id}) {
		foreach my $y (@{$pass_paf{$id}}) {
			print $F_PAF_OUTPUT "$y\n";
		}
	}
	close $F_PAF_OUTPUT;
	
	my (@sam_header, %sam_reads);
	open my $F_SAM_INPUT, '<', $sam or die "### $prog_id ERROR: can't open sample SAM for sample $s_name!\n";
	while (<$F_SAM_INPUT>) {
		chomp;
		if ($_ =~ m/^\@/) {
			push @sam_header, $_;
		} else {
			my $sam_body = $_;
			my @s = split '\t', $sam_body;
			my $sam_id = $s[0];
			if (exists $pass_id{$sam_id}) {
				push @{$sam_reads{$sam_id}}, $sam_body;
			}
		}
	}
	close $F_SAM_INPUT;
	
	open my $F_SAM_OUTPUT, '>', $out_sam or die "### $prog_id ERROR: an't write filtered SAM for sample $s_name!\n";
	print $F_SAM_OUTPUT "$_\n" for @sam_header;
	foreach my $id (@pass_id_sorted) {
		foreach my $r (@{$sam_reads{$id}}) {
			print $F_SAM_OUTPUT "$r\n";
		}
	}
	close $F_SAM_OUTPUT;
}

my ($f_io_summary, $q_threshold);
if (@ARGV) {
	($f_io_summary, $q_threshold) = @ARGV;
	die "### $prog_id ERROR: I/O summary cannot be located:\n\t$f_io_summary" unless (-e $f_io_summary);
	unless (defined $q_threshold) {
		printf STDERR "### $prog_id WARNING: threshold not defined, assumes 50.\n";
		$q_threshold = 50; 
	}
} else {
	die "quality_filter_paf_sam <IO summary file> <numeric quality threshold>\n";
}

my (%header, %samples, %paf_hq, %paf_lq, %sam_hq, %sam_lq, %paf_out_hq, %paf_out_lq, %sam_out_hq, %sam_out_lq);
open my $F_SUMMARY, '<', $f_io_summary or die "### $prog_id ERROR: can't read I/O summary file!\n";
my $n_cur = 0;
while (<$F_SUMMARY>) {
	chomp;
	my @s = split ',', $_;
	if ($n_cur == 0) {
		%header = %{get_index(\@s, 
			qw(minimap2_hq_paf minimap2_lq_paf minimap2_hq_sam minimap2_lq_sam qfilt_hq_paf qfilt_lq_paf qfilt_hq_sam qfilt_lq_sam))};
	} else {
		my $id = $s[0];
		$samples{$id}++;
		$paf_hq{$id}     = $s[$header{minimap2_hq_paf}] if (defined $header{minimap2_hq_paf}) ;
		$paf_lq{$id}     = $s[$header{minimap2_lq_paf}] if (defined $header{minimap2_lq_paf}) ;
		$sam_hq{$id}     = $s[$header{minimap2_hq_sam}] if (defined $header{minimap2_hq_sam}) ;
		$sam_lq{$id}     = $s[$header{minimap2_lq_sam}] if (defined $header{minimap2_lq_sam}) ;
		$paf_out_hq{$id} = $s[$header{qfilt_hq_paf}] if (defined $header{qfilt_hq_paf}) ;
		$paf_out_lq{$id} = $s[$header{qfilt_lq_paf}] if (defined $header{qfilt_lq_paf}) ;
		$sam_out_hq{$id} = $s[$header{qfilt_hq_sam}] if (defined $header{qfilt_hq_sam}) ;
		$sam_out_lq{$id} = $s[$header{qfilt_lq_sam}] if (defined $header{qfilt_lq_sam}) ;
	}
	$n_cur++;
}
close $F_SUMMARY;

foreach my $id (keys %samples) {
	printf STDERR "### $prog_id: applying quality filters to HQ seq data for sample %s...", $id;
	qfilt($paf_hq{$id} , $sam_hq{$id}, $paf_out_hq{$id}, $sam_out_hq{$id}, $q_threshold, $id);
	print STDERR "\n";
	if ((exists $paf_lq{$id}) && (exists $sam_lq{$id}) && (exists $paf_out_lq{$id}) && (exists $sam_out_lq{$id})) {
		printf STDERR "### $prog_id: applying quality filters to LQ seq data for sample %s...", $id;
		qfilt($paf_lq{$id} , $sam_lq{$id}, $paf_out_lq{$id}, $sam_out_lq{$id}, $q_threshold, $id);
		print STDERR "\n";
	}
}

exit 0;
