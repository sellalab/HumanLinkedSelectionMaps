#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <time.h>

#include "config.h"
#include "util.h"
#include "dist.h"
#include "rectab.h"

//#include <genome_db.h>

//#include <analysis.h>

#include "bkgd.h"
#include "bkgd_param.h"
#include "bkgd_interp.h"
#include "bkgd_files.h"



/**
 * Retrieves a table of recombination rates for current chromosome
 * from database
 */
static RecRateTable *get_rectab(Config *config, Chromosome *chr) {
  long n_sf;
  //  SeqCoord *region;
  //  SeqFeatureDBA *sf_dba;
  SeqFeature *sf;
  char *rec_tab_name;
  double rec_scale;
  RecRateTable *rtab;
  
  //  region = analysis_get_region(an);
  //  sf_dba = seqfeat_dba_create(an->dbc);

  rec_tab_name = config_get_str(config, "RECOMB_RATE_TABLE");
  rec_scale = config_get_double(config, "RECOMB_RATE_SCALE");

  /* get recombination rates from database */
  sf = load_genetic_map_from_file(rec_tab_name, chr, &n_sf);
  //  sf = seqfeat_dba_fetch_by_region(sf_dba, region, rec_tab_name, &n_sf);
  rtab = rectab_from_feats(sf, n_sf, chr->len, rec_scale);
  seqfeat_array_free(sf, n_sf);

  return rtab;
}





/**
 * Converts provided recombination rates to recombination distances
 * along chromosome and retrieves a list of rec-dist coords
 * representing conserved blocks. Conserved blocks are split when
 * recombination rates change in order to make integral approximation
 * of background selection sum more efficient.
 */
static GList *get_cons_rec_dists(Config *config, RecRateTable *rtab) {
  long i, j, n_sf, ttl_left, ttl_right, ttl_cons, last_end;
  //  long len; 
  //  SeqCoord *region;
  //  SeqFeatureDBA *sf_dba;
  SeqFeature *sf;
  ConsBlock *cblk;
  GList *cons_list, *cur;
  char *cons_tab;
  
  //  region = analysis_get_region(an);
  //  len = region->end - region->start + 1;
  //  sf_dba = seqfeat_dba_create(an->dbc);
  cons_tab = config_get_str(config, "CONS_TABLE");

  /* get conserved elements from database */
  sf = load_conserved_from_file(cons_tab, config_get_str(config, "CHROMOSOME_NAME"), &n_sf);
  //  sf = seqfeat_dba_fetch_by_region(sf_dba, region, cons_tab, &n_sf);
  //  /* order conserved elements */
  //  qsort(sf, n_sf, sizeof(SeqFeature), seqfeat_cmp_nostrand);

  /* Create list of conserved block coords splitting elements when
   * they span a change in recombination rate. The coordinates are
   * initially set to physical positions, and then updated to be
   * recombination distances.
   */
  cons_list = NULL;
  ttl_cons = 0;
  last_end = -1;

  for(i = 0; i < n_sf; i++) {
    cblk = g_new(ConsBlock, 1);
    cblk->start = sf[i].c.start;

    /* sanity check: we don't want overlapping blocks */
    // DEBUG: allow overlap where new_start = old_end only!!
    // if(cblk->start <= last_end) {
    if(cblk->start < last_end) {
      g_error("get_cons_rec_dists: conserved blocks should not overlap: "
	      "cur_start=%ld, last_end=%ld", cblk->start, last_end);
    }
    
    // add length of each segment to total conservered tally
    ttl_cons += sf[i].c.end - sf[i].c.start + 1;
    cblk->r = rectab_rate(rtab, sf[i].c.start);

    /* re-write conserved blocks in genetic map coordinates
     * iterate over the entire length of each block and
     * split the block if there is a change in genetic
     * map coordinates at a given point
     */ 
    for(j = sf[i].c.start; j < sf[i].c.end; j++) {
      if(rectab_rate(rtab,j+1) != cblk->r) {   	
        /* change in rate, end previous block */
      	cblk->end = j;
      	last_end = cblk->end;

        // DEBUG: write the new truncated segment to debug file
        // fprintf(db_blk_fh, "%ld %ld\n", cblk->start, cblk->end);

        // append previous block to running list
      	cons_list = g_list_append(cons_list, cblk);
      	
        /* 	fprintf(stderr, "splitting cons segment %ld-%ld into " */
        /* 		"%ld-%ld and %ld-%ld (rec-rate change from %g to %g)\n", */
        /* 		cblk->start, sf[i].c.end, */
        /* 		cblk->start, cblk->end, */
        /* 		j+1, sf[i].c.end, cblk->r, rectab_rate(rtab,j+1)); */

      	/* start new block */
      	cblk = g_new(ConsBlock, 1);
      	cblk->start = j+1;
      	cblk->r = rectab_rate(rtab, j+1);
      }
    }
    /* end current block */
    cblk->end = sf[i].c.end;
    last_end = cblk->end;
    cons_list = g_list_append(cons_list, cblk);
  }

  // DEBUG: print total number of conserved bases stderr
  fprintf(stderr, "total conserved sites summed over conserved regions: %ld\n", ttl_cons);
  seqfeat_array_free(sf, n_sf);
  
  /* convert cons block phys coords to rec-dist coords */
  ttl_left = 0;
  ttl_right = ttl_cons;

  // DEBUG: debug file path
  // char *chname, *db_rblk_dir, *db_rblk_path;
  // chname = config_get_str(config, "CHROMOSOME_NAME");
  // db_rblk_dir = config_get_str(config, "OUTPUT_DIR");
  // db_rblk_path = g_strconcat(db_rblk_dir, chname, ".rstart.rend.bed", NULL);
  // FILE *db_rblk_fh;
  // db_rblk_fh = fopen(db_rblk_path, "w");

  // DEBUG: cons_list2 is for the second split on genetic length of segments
  // GList *cons_list2;
  // cons_list2 = NULL;
  // ConsBlock *sub_cblk;
  // long slen, nseg, blk_start, blk_end, blk_len;
  // double blk_r_start, blk_r_end, blk_rlen, tot_segs;

  cur = cons_list;
  while(cur != NULL) {
    cblk = cur->data;
    cblk->r_start = rectab_rpos(rtab, cblk->start);
    cblk->r_end = rectab_rpos(rtab, cblk->end);
    /* keep track of number of cons sites to left and right of each
     * blk (not including the block itself)
     */
    cblk->left_ttl = ttl_left;
    ttl_left  += cblk->end - cblk->start + 1;
    ttl_right -= cblk->end - cblk->start + 1;
    cblk->right_ttl = ttl_right;
   
    // DEBUG: get the physical and genetic length of the block
    // blk_start = cblk->start;
    // blk_end = cblk->end;
    // blk_len = blk_end - blk_start;
    // blk_r_start = cblk->r_start;
    // blk_r_end = cblk->r_end;
    // blk_rlen = blk_r_end - blk_r_start;

    // DEBUG: for blocks greater than the average genetic length, split into sub-blocks
    // if (FALSE) {
    //   // if (blk_rlen > 1e-08) {
    //   /* num sub-blocks = min(ceil(blk_rlen / 1e-08), blk_len) */
    //   if (ceil(blk_rlen / 1e-08) < blk_len)
    //     nseg = ceil(blk_rlen / 1e-08);
    //   else
    //     nseg = blk_len;
    //   /* find max size for blocks that will be <= 1e-08 (except when limit=1bp is reached) */
    //   slen = blk_len / nseg;
    //   /* true number of segments + remainder with new segment size */
    //   nseg = ceil((double)(blk_len) / (double)(slen));

    //   // DEBUG print
    //   if ((nseg*slen - blk_len) >= slen) {
    //     fprintf(db_rblk_fh, "# ERROR: ");
    //     fprintf(db_rblk_fh, "SPLIT: nseg=%ld. slen=%ld. tot=%ld. init=%ld. remain=%ld\n", nseg, slen, nseg*slen, blk_len, (nseg*slen)-blk_len);
    //     fflush(db_rblk_fh);
    //     }
    //   else {
    //     fprintf(db_rblk_fh, "# SPLIT: nseg=%ld. slen=%ld. tot=%ld. init=%ld. remain=%ld\n", nseg, slen, nseg*slen, blk_len, (nseg*slen)-blk_len);
    //     fflush(db_rblk_fh);
    //   }

    //   /* create the first n-1 sub-blocks; the final block may need to be truncated */
    //   for (i=1; i<nseg; i++) {
    //     tot_segs+=1;
    //     sub_cblk = g_new(ConsBlock, 1);
    //     /* set start-end coords */
    //     sub_cblk->start = blk_start;
    //     sub_cblk->end = blk_start + slen;
    //     /* set rate to the current block (uniform for all sub-blocks) */
    //     sub_cblk->r = cblk->r;
    //     /* set the genetic start and end points to sub-block */
    //     sub_cblk->r_start = rectab_rpos(rtab, sub_cblk->start);
    //     sub_cblk->r_end = rectab_rpos(rtab, sub_cblk->end);
    //     /* update conserved counts to right and left */
    //     sub_cblk->left_ttl = ttl_left;
    //     ttl_left += sub_cblk->end - sub_cblk->start + 1;
    //     ttl_right -= sub_cblk->end - sub_cblk->start + 1;
    //     sub_cblk->right_ttl = ttl_right;
    //     /* append sub-block */
    //     cons_list2 = g_list_append(cons_list2, sub_cblk);
    //     /* update next sub-block start point */
    //     blk_start = sub_cblk->end;

    //     // DEBUG print
    //     fprintf(db_rblk_fh, "%ld\t%ld\t%.6g\t%.6g\tsplit_%ld\n", sub_cblk->start, sub_cblk->end, sub_cblk->r_start, sub_cblk->r_end, i);
    //     fflush(db_rblk_fh);
    //   }
    //   assert (blk_start < blk_end);
    //   /* compete the final sub-block and truncate length if necessary */
    //   tot_segs+=1.0;
    //   sub_cblk = g_new(ConsBlock, 1);
    //   /* set start-end coords */
    //   sub_cblk->start = blk_start;
    //   sub_cblk->end = blk_end;
    //   /* set rate to the current block (uniform for all sub-blocks) */
    //   sub_cblk->r = cblk->r;
    //   /* set the genetic start and end points to sub-block */
    //   sub_cblk->r_start = rectab_rpos(rtab, sub_cblk->start);
    //   sub_cblk->r_end = rectab_rpos(rtab, sub_cblk->end);
    //   /* update conserved counts to right and left */
    //   sub_cblk->left_ttl = ttl_left;
    //   ttl_left += sub_cblk->end - sub_cblk->start + 1;
    //   ttl_right -= sub_cblk->end - sub_cblk->start + 1;
    //   sub_cblk->right_ttl = ttl_right;
    //   /* append sub-block */
    //   cons_list2 = g_list_append(cons_list2, sub_cblk);

    //   // DEBUG print
    //   fprintf(db_rblk_fh, "%ld\t%ld\t%.6g\t%.6g\tsplit_end_%ld\n", sub_cblk->start, sub_cblk->end, sub_cblk->r_start, sub_cblk->r_end, i);
    //   fflush(db_rblk_fh);
    // }

    // else {
    //   // tot_segs+=1.0;
    //   /* keep track of number of cons sites to left and right of each
    //    * blk (not including the block itself)
    //    */
    //   cblk->left_ttl = ttl_left;
    //   ttl_left  += cblk->end - cblk->start + 1;
    //   ttl_right -= cblk->end - cblk->start + 1;
    //   cblk->right_ttl = ttl_right;
    //   /* append sub-block */
    //   cons_list2 = g_list_append(cons_list2, cblk);

    //   // DEBUG printout
    //   // fprintf(db_rblk_fh, "%ld\t%ld\t%.6g\t%.6g\tpassed\n", cblk->start, cblk->end, cblk->r_start, cblk->r_end);
    //   // fflush(db_rblk_fh);
    // }
  cur = g_list_next(cur);
  }

  // fprintf(stderr, "finished splitting conserved segments into %g segments\n", tot_segs);
  // fflush(stderr);
  // fclose(db_rblk_fh);  // DEBUG: close rblk log file

  return cons_list;
  // DEBUG: return the modifed cons list with split blocks
  // return cons_list2;
}



/* 
 * Calculates background selection strength at each position along a
 * chromosome using provided parameters, recombination map and list of
 * conserved elements.
 * 
 */
// DEBUG: add a new param for debug file handle in this definition
// void calc_bkgd_chr(Config *config, Chromosome *chr, RecRateTable *rtab, GList *cons_list,
            // BkgdParam *parm,  FILE *out_fh, FILE *debug_fh) {

void calc_bkgd_chr(Config *config, Chromosome *chr, RecRateTable *rtab, GList *cons_list,
            BkgdParam *parm,  FILE *out_fh) {

  double b;
  long pos, b_len;
  int b_int, prev_b_int;
  GList *next_cons;
  BkgdInterp *bgi;
  next_cons = cons_list;

  /* b is background selection strength */
  prev_b_int = b_int = -1;

  // length of segment where b within +/- epsilon
  b_len = 0;

  /* Create interpolator to estimate B values at positions along chr */
  // bgi = bkgd_interp_new(rtab, chr->len, cons_list, parm);
  // DAVID: init bgi with chr
  bgi = bkgd_interp_new(rtab, chr, cons_list, parm);

  // pos = 1;
  // DAVID: begin at the first position in the partition
  pos = chr->p_start;

  // DEBUG:
  fprintf(stderr, "initial pos=%ld\n", pos);
  fflush(stderr);

  // DAVID: 1% of partition length
  long pct;
  pct = (chr->p_end - chr->p_start + 1) / 100;


  // DEBUG:
  long n_segs, seg_sum;
  n_segs = 0;
  seg_sum = 0;

  // while(pos <= chr->len) {  
  // DAVID: only calculate b until the end of the partition
  while(pos <= chr->p_end) {
    b = bkgd_interp_eval(bgi, pos);
    // DEBUG: pass debug file handle to interpolation function to record points, bvals
    // b = bkgd_interp_eval(bgi, pos, debug_fh);
    
    // if((pos % 1000000)==0) {
    // DAVID: count 1% increments of partition length
    if ((pos % pct) == 0) {
      fprintf(stderr, ".");
      fflush(stderr);
    }
    /* truncate to integer for range determined by BKGD_SCALE */
    // b_int = (int)floor(b * BKGD_SCALE + 0.5);

    /* DAVID: change here needed for printing little b's on a new grid
     * For the little b regime scale DIVIDING b by some scale
     * factor -- e.g., if we want 100 bins, we divide by the 
     * bin size and round.
     */
    double eps_b, eps_B;
    eps_B = 1.0 / (double)BKGD_SCALE;
    eps_b = -gsl_log1p(-eps_B);
    b_int = (int)floor(-b / eps_b + 0.5);

    /* only print out value if rounded value is different from prev one */
    if(prev_b_int != b_int) {
      if(prev_b_int >= 0) {
        fprintf(out_fh, "%d %ld\n", prev_b_int, b_len);
      }
      // DEBUG:
      n_segs += 1;
      seg_sum += b_len;

      prev_b_int = b_int;
      b_len = 0;
    }
    pos++;
    b_len++;
  }

  /* print out final value */
  long plen;  // get partition length for checking final length
  int diff; 
  plen = chr->p_end - chr->p_start + 1;
  if(b_len > 0) {
    // DEBUG:
    n_segs += 1;
    seg_sum += b_len;
    // get the difference between segment sum and partition length (should be 0)
    diff = plen - seg_sum;
    // add the difference to blen to close any gaps
    b_len += diff;
    fprintf(out_fh, "%d %ld\n", b_int, b_len);
    /* fprintf(stderr, "%ld %d %d\n", pos, prev_b_int, b_len); */
  }

  fprintf(stderr, "\n");

  // DEBUG:
  fprintf(stderr, "n_segs=%ld, seg_sum=%ld, pos=%ld, b_len=%ld\n", n_segs, seg_sum, pos, b_len);
  fflush(stderr);

  bkgd_interp_free(bgi);

}



int main(int argc, char **argv) {
  Config *config;
  Chromosome *chr;
  GList *cons_list, *cur;
  char *out_dir, *out_token;
  // DEBUG: declare debug file path
  // char *out_dir, *fout, *out_token, *debug_path;
  RecRateTable *rtab;
  // DEBUG: declare debug file handle
  // FILE *out_fh, *debug_fh;
  FILE *out_fh;
  BkgdParam *parm;
  
  if(argc != 3) {
    fprintf(stderr, "usage: %s <config_file> <output_file_name>\n",
	    argv[0]);
    exit(2);
  }

  // record the time when the program begins
  time_t time_0 = time(NULL);

  fprintf(stderr, "Reading config\n");
  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);

  /* GUY: Changed the length to be a direct configuration input
   *  chr_features_file = config_get_str(config, "CHROMOSOME_FEATURES");
   *  get_chr_features(chr, argv[2], chr_features_file);
   */
  chr = g_new(Chromosome, 1);
  chr->len = atol(config_get_str(config, "CHROMOSOME_LENGTH"));
  chr->id = -1;
  chr->assembly = NULL;
  chr->name = g_new(char, 256);
  sscanf(config_get_str(config, "CHROMOSOME_NAME"), "%s", chr->name);

  /* DAVID: Partition chromosome into regions for parallelizing 
   * large chroms over multiple processes. Use length and index
   * of partition to get region (start, end).
   */
  long pidx, plen;
  pidx = atol(config_get_str(config, "CHR_PARTITION_INDEX"));
  plen = (long)atof(config_get_str(config, "CHR_PARTITION_LENGTH"));
  fprintf(stderr, "CHR_PARTITION_INDEX=%ld\nCHR_PARTITION_LENGTH=%ld\n", pidx, plen);
  fflush(stderr);
  // bypass partitioning option if pidx < 0 or plen <= 0
  if ((pidx < 0) || (plen == 0)) {
    chr->p_start = 1;
    chr->p_end = chr->len;
  }
  // otherwise set partition coordinates
  else {
    // set region start; must be within the current chromosome
    chr->p_start = pidx * plen + 1;
    if (chr->p_start >= chr->len) {
      fprintf(stderr, "ERROR: p_start=%ld exceeds %s length\n", chr->p_start, chr->name);
      exit(2);
    }
    // set region end; bound by chrom length
    chr->p_end = long_min(chr->len, chr->p_start + plen - 1);
  }
  fprintf(stderr, "CHR_PARTITION_START=%ld\nCHR_PARTITION_END=%ld\n", chr->p_start, chr->p_end);
  fflush(stderr);

  out_dir = config_get_str(config, "OUTPUT_DIR");
  // GUY: changed filneame to be a commandline argument
  // out_token = config_get_str(config, "OUTPUT_TOKEN");
  out_token = argv[2];

  // new parameter struct
  parm = bkgd_param_new(config);

  // open bmap file
  char fout[1000];
  snprintf(fout, 999, "%s%s.bkgd", out_dir, out_token);
  out_fh = fopen(fout, "w");

  // DEBUG LINES: create a new file for writing interpolation box points and exact B values
  // debug_path = g_strconcat(out_dir, out_token, ".exactB.txt", NULL);
  // debug_fh = fopen(debug_path, "w");

  if(out_fh == NULL) {
    g_error("main: could not open output file '%s'", fout); 
  }
  // g_free(fout);

  // write input params to header of bmap file
	print_params(out_fh, chr, argv[1]);

  // new recombination map struct
  fprintf(stderr, "retrieving recombination rates\n");
  rtab = get_rectab(config, chr);

  // get recombination distances of conserved sites, and convert rec rates to dists
  fprintf(stderr, "calculating conserved rec dists\n");
  cons_list = get_cons_rec_dists(config, rtab);
  fprintf(stderr, "total recomb dist for %s: %gM\n", chr->name, rtab->chr_r_len);

  // record the time to load annotations and build lookup tables. reset time_0
  double time_startup = difftime(time(NULL), time_0) / 60.0;
  fprintf(stderr, "startup time (minutes) = %.3f\n", time_startup);
  time_0 = time(NULL);

  // calculate strength of background selection at each position on chr
  fprintf(stderr, "calculating b vals\n");
  calc_bkgd_chr(config, chr, rtab, cons_list, parm, out_fh);
  // DEBUG: pass the debug file handle to calc_bkgd_chr
  // calc_bkgd_chr(config, chr, rtab, cons_list, parm, out_fh, debug_fh);

  // record the time spend calculating b values across the chromosome/partition. reset time_0
  double time_bmap = difftime(time(NULL), time_0) / 60.0;
  fprintf(stderr, "b calculation time (minutes) = %.3f\n", time_bmap);
  time_0 = time(NULL);

  // free various malloc blocks, close files
  fprintf(stderr, "freeing recombination distances\n");
  rectab_free(rtab);    

  fprintf(stderr, "freeing conserved blocks\n");
  cur = cons_list;
  while(cur != NULL) {
    g_free(cur->data);
    cur = g_list_next(cur);
  }
  g_list_free(cons_list);

  fclose(out_fh);
  // DEBUG: close the debug file
  // fclose(debug_fh);
  
  fprintf(stderr, "freeing parameters and interpolation tables\n");
  bkgd_param_free(parm);

  fprintf(stderr, "freeing analysis and config\n");
  config_free(config);
  
  // record the shutdown time after after bmap has been completed
  double time_shutdown = difftime(time(NULL), time_0) / 60.0;
  fprintf(stderr, "shutdown time (minutes) = %.3f\n", time_shutdown);

  return 0;  
}

