#include <glib.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <assert.h>

#include "config.h"
#include "numer.h"
#include "bkgd_interp.h"
#include "bkgd_param.h"
#include "bkgd_point.h"
#include "bkgd.h"
#include "bkgd_intg.h"


/* maximum distance between interpolation points in morgans */
double BKGD_INTERP_MAX_DIST = 0.0001;
//#define BKGD_INTERP_MAX_DIST 1e-4

/* if the bkgd parameter changes by this much between chosen
 * interpolation points, find closer points to use
 */
double BKGD_CHANGE_THRESH = 0.002;
//#define BKGD_CHANGE_THRESH 0.002

/* when doing quadratic extrapolation aim for this much change in B */
double BKGD_OFFSET = 0.001;
//#define BKGD_OFFSET 0.001





/**
 * Returns the location of the next sample that we want to take given
 * the current location by fitting a quadratic at the current
 * position.
 */
static double get_r_next(double x, double bkgd, double b_drv1, double b_drv2) {
  double a, b, c, c_new, discr, sq, x1,x2,x3,x4, closest;

  /* Fit a curve to B at x by finding a with the same value, 1st and
   * 2nd derivates at this point.
   */

  /* find coefficients of quadtratic */
  a = 0.5 * b_drv2;
  b = b_drv1 - b_drv2 * x;
  c = a*x*x - b_drv1*x + bkgd;


  /* now we want to find point where we estimate B will have changed
   * by BKGD_OFFSET offset using the quadratic equation
   * DAVID: this offset will have to be set in little B log-space
   */

  /* find next point that is offset GREATER */
  c_new = c - bkgd-BKGD_OFFSET;
  discr = b*b - 4*a*c_new;

  if(discr < 0.0) {
    x1 = x2 = x + BKGD_INTERP_MAX_DIST;
  } 
  else {
    sq = sqrt(discr);
    x1 = (-b + sq)/(2.0*a);
    x2 = (-b - sq)/(2.0*a);
  }
  
  /* find next point that is offset LOWER */
  c_new = c - bkgd+BKGD_OFFSET;

  discr = b*b - 4*a*c_new;

  if(discr < 0.0) {
    x3 = x4 = x + BKGD_INTERP_MAX_DIST;
  } 
  else {
    sq = sqrt(discr);
    x3 = (-b + sq)/(2.0*a);
    x4 = (-b - sq)/(2.0*a);
    /* FIXED: 
     * changed x4 = (-b + sq)/(2.0*a); to x4 = (-b - sq)/(2.0*a);
     * (David)
     */
  }

  /* take closest point which is greater than x */
  closest = x+BKGD_INTERP_MAX_DIST;
  if(x1 > x && x1 < closest) {
    closest = x1;
  }
  if(x2 > x && x2 < closest) {
    closest = x2;
  }
  if(x3 > x && x3 < closest) {
    closest = x3;
  }
  if(x4 > x && x4 < closest) {
    closest = x4;
  }

  
  /*   fprintf(stderr, "x1=%g, x2=%g, x3=%g, x4=%g closest=%g\n",  */
  /* 	  x1, x2, x3, x4, closest); */

  return closest;
}


/**
 * Decides on what position should be evaluated next using quadratic
 * extrapolation from last_p, and sets the attributes of p
 * appropriately.
 */
static void get_next_point(BkgdInterp *bgi, BkgdPoint *last_p, BkgdPoint *p) {
  double next_r, r_delta, b_delta;
  long pos_delta;
  int keep_point;
  ConsBlock *cblk;
  BkgdPoint *saved_p;
  double next_r_cons, next_r_mid;


  /* predict a good next r position */
  next_r = get_r_next(last_p->r_pos,  last_p->b, last_p->b_drv1, last_p->b_drv2);
  r_delta = 0.0;


  // BEGIN MODIFIED POINT SEARCH
  // ***************************

  /*
   * DESCRIPTION:
   * use the site 1bp left of next cons to bound step size between points
   * interpolate the genetic map position of the cons-1 site and set next_r
   * to be the minimum of the genetic map position and the get_r_next() output

   * if the last_p is cons-1, automatically jump to the next neutral
   * site 1bp right of next cons -- this may require advancing next cons
   * through contiguous conserved blocks for blocks that were split
   * where recombination rate changed within the block

   * try current cons block end+1 for the next neutral site
   * if one or more of the NEXT cons blocks are contiguously
   * joined as a result of genetic map splitting, end+1 will 
   * overlap the start of the next conserved block. in this
   * case, we advance the next conserved block again and try
   * end+1 for the next conserved site -- the process is repeated
   * until we've jumped over all contigous cons blocks to reach
   * the next neutral site

   * jump through regions of no recombination by advancing next cons until
   * the genetic position of cblk->start increases and then take one step 
   * back leaving 2 points of identical b values and r_pos bracketing each
   * region of no recombination
   *
   * NOTE: this whole section should be moved to a separate function
   */

  // use enforced point stops if flagged
  if (bgi->parm->enforce_bgi_pts) {

    // check that there are conserved sites ahead, otherwise use standard algorithm
    if (bgi->next_cons != NULL) {

      // get next conserved block associated with last_p
      cblk = bgi->next_cons->data;

      // jump automatically if last_p is at cons-1 or it falls inside a cons block
      // NOTE: technically we should only fall inside a cons block if position 1 in chrom is in cons block
      if ((last_p->pos >= cblk->start-1) && (last_p->pos <= cblk->end)) {

        /* JUMPING TO NEXT NEUTRAL SITE:
         * try the position 1bp right of the end of the current block
         * advance to the next conserved block
         * break if there are no more blocks
         * if there are more blocks, update cblk with new block
         * break if new block start > last block end+1
         */
        while (TRUE) {
          p->pos = cblk->end+1;
          bgi->next_cons = g_list_next(bgi->next_cons);
          if (bgi->next_cons == NULL)
            break;
          cblk = bgi->next_cons->data;
          if (cblk->start > p->pos)
            break;
        }
        // store this value in the interpolator for calculating cblk midpoints later
        // NOTE: even if we backup to stay in the partition, this value should still be w.r.t. cons blocks
        // bgi->last_cblk_rend = rectab_rpos(bgi->rtab, p->pos);

        // EXCEPTIONAL: if final cblk extends through last bp of the chrom, use last chrom position
        if (p->pos > bgi->chr_len)
          p->pos = bgi->chr_len;
        // PARTITIONED CHROM VERSION: if pos jumped out of partition, backup to last site of partition
        // if (p->pos > bgi->chr->p_end)
        //   p->pos = bgi->chr->p_end;

        // get genetic map position of the final point
        p->r_pos = rectab_rpos(bgi->rtab, p->pos);
        bgi->last_cblk_rend = p->r_pos;

        // calculate exact b at the new point
        bkgd_calc_b(p, bgi->cons_list, bgi->next_cons, bgi->parm, bgi->rtab->chr_r_len);

        // calculate the new slope
        r_delta = (p->r_pos - last_p->r_pos);
        if(r_delta == 0.0)
          bgi->slope = 0.0;
        else
          bgi->slope = (p->b - last_p->b) / r_delta;

        // bgi->return_condition = "jump over group of contiguous cons blocks";
        // bgi->return_condition = "R1";
        return;
      }

      // if the first condition isn't triggered than cblk is fully right of pos
      assert (last_p->pos < cblk->start-1);

      // get the genetic map position of cons block to right start-1
      next_r_cons = rectab_rpos(bgi->rtab, cblk->start-1);

      // this MUST be nonnegative
      assert (next_r_cons - bgi->last_cblk_rend >= 0);

      // get the genetic map midpoint between cons block to left end+1 and cons block to right start-1
      next_r_mid = bgi->last_cblk_rend + ((next_r_cons - bgi->last_cblk_rend) / 2.0);

      // if the last point was behind the midpoint, use midpoint as max for next_r
      if ((last_p->r_pos < next_r_mid) && (next_r_mid < next_r))
        next_r = next_r_mid;

      // if last point was >= the midpoint, use cons block start-1 as max for next_r
      if (next_r > next_r_cons)
        next_r = next_r_cons;

      // if the genetic coords haven't changed then the region has no recombination and we jump
      // over conserved segments until recombination resumes
      if (next_r == last_p->r_pos) {

        // we assume the new position is at min 1bp above previous position
        p->pos = last_p->pos+1;

        // increment p->pos by 1bp until p->r_pos > last_p->r_pos, then backup 1 so we arrive at the furthest
        // physical position from last_p->pos where the map distance stays the same
        while (TRUE) {
          p->pos += 1;
          // if end of chrom reached, break
          if (p->pos > bgi->chr_len) {
            p->pos -= 1;
            break;
          }
          // update next cons if we pass beyond the current cblk end+1
          // (+1 buffers against the subtraction at breakpoint)
          if (p->pos > cblk->end+1) {
            bgi->next_cons = g_list_next(bgi->next_cons);
            if (bgi->next_cons != NULL)
              cblk = bgi->next_cons->data;
          }
          // as soon as rpos changes backup 1 and break
          p->r_pos = rectab_rpos(bgi->rtab, p->pos);
          if (p->r_pos > last_p->r_pos) {
            p->pos -= 1;
            break;
          }
        }

        // final check if we are beyond current cblk after adjusting back -1
        if (p->pos > cblk->end) {
          bgi->next_cons = g_list_next(bgi->next_cons);
          if (bgi->next_cons != NULL)
            cblk = bgi->next_cons->data;
        }

        // jump to the first cons block where next_r increases if there are any left in the list
        // while (TRUE) {
        //   bgi->next_cons = g_list_next(bgi->next_cons);
        //   if (bgi->next_cons == NULL)
        //     bgi->next_cons = g_list_last(bgi->cons_list);
        //     break;
        //   cblk = bgi->next_cons->data;
        //   next_r = rectab_rpos(bgi->rtab, cblk->start-1);
        //   if (next_r > last_p->r_pos)
        //     bgi->next_cons = g_list_previous(bgi->next_cons);
        //     break;
        // }

        // reset current cons block
        // cblk = bgi->next_cons->data;

        // update physical position to minimum of p_end vs. cblk-1
        // if (cblk->start-1 > bgi->chr->p_end)
        //   p->pos = bgi->chr->p_end;
        // else
        //   p->pos = cblk->start-1;


        // genetic position has not changed so neither has b or its derivatives
        p->r_pos = last_p->r_pos;
        p->b = last_p->b;
        p->b_drv1 = last_p->b_drv1;
        p->b_drv2 = last_p->b_drv2;

        // slope is zero since there was no change to the genetic map position
        bgi->slope = 0.0;

        // bgi->return_condition = "jump over conserved blocks in a 0 M/bp region";
        // bgi->return_condition = "R2";
        return;
      }
    }
  }
  // END MODIFIED POINT SEARCH
  // *************************


  /* check the queue of saved positions before evaluating B at new positions */
  if(bgi->p_queue->length > 0) {
    saved_p = g_queue_peek_head(bgi->p_queue);

    if(saved_p->r_pos <= next_r) {
      /* saved position is closer than predicted position, so try
       * it instead
       */

      /* calculate differences in saved position and last position */
      b_delta   = fabs(saved_p->b - last_p->b);
      pos_delta = saved_p->pos - last_p->pos;
      r_delta   = (saved_p->r_pos - last_p->r_pos);

      if(pos_delta > 1 && b_delta > BKGD_CHANGE_THRESH && r_delta > 0.0) {
      	/* Don't use saved position because still too far away */
      	/* instead try a position that is 1/4 as far away */
      	// next_r = last_p->r_pos + r_delta * 0.25;
        // DEBUG: try taking less conservative step 1/2 r_delta
        next_r = last_p->r_pos + r_delta * 0.5;
      	keep_point = FALSE;
      } 
      else {
      	/* saved position is good */
      	g_queue_pop_head(bgi->p_queue);

        /*  	fprintf(stderr, "using saved point at %g, queue len=%d\n", */
        /*  		saved_p->r_pos, bgi->p_queue->length); */
        /* 	fprintf(stderr, "last_p->b=%g, saved_p->b: %g, " */
        /* 		"b_delta=%g, pos_delta=%ld, r_delta=%g\n", */
        /* 		last_p->b, saved_p->b, b_delta, pos_delta, r_delta); */

      	bkgd_point_copy(saved_p, p);
      	g_free(saved_p);
      	keep_point = TRUE;
      }
    } 
    else {
      /* predicted position is closer, try it instead */
      keep_point = FALSE;
    }
  } 
  else {
    /* there are no saved positions so try predicted position */
    keep_point = FALSE;
  }

  /* keep trying closer points until there are no closer points,
   * or the change threshold is within our tolerance
   */
  while(!keep_point) {
    p->pos = last_p->pos+1;

    /* find first base that is greater than next r genetic map position */
    while((p->pos < bgi->chr_len) && (rectab_rpos(bgi->rtab, p->pos) < next_r)) {
    // DAVID: bound p->pos by partition end
    // while((p->pos < bgi->chr->p_end) && (rectab_rpos(bgi->rtab, p->pos) < next_r)) {
      p->pos += 1;
    }

    if((p->pos > last_p->pos+1) && (rectab_rpos(bgi->rtab,p->pos) > next_r)) {
      /* backup one base so we don't exceed r pos we were aiming for */
      p->pos -= 1;
    }
  
    p->r_pos = rectab_rpos(bgi->rtab,p->pos);


    /* 
     * NOTE ON THE MODIFIED VERSION: in theory this step should never 
     * be triggered because we only move one conserved block at a time
     */

    /* Update next_cons so that it points to next conserved block.
     * Remember that p can move backwards if we chose a site that
     * was too far ahead. When this happens we may have to backup in
     * the conserved block list instead of moving forwards.
     */

    /* first backup */
    while(bgi->next_cons != NULL) {
      cblk = bgi->next_cons->data;
      // NOTE: since we step up to cblk->start-1 at MOST, 
      // this condition should (never?) be true if we
      // advance next_cons inside the cblk->end+1 condition
      if(cblk->start <= p->pos) {
	      break;
      }
      bgi->next_cons = g_list_previous(bgi->next_cons);
    }

    /* if next_cons was pointing at the first cblk, the
    first backup may have moved its pointer to NULL */
    if(bgi->next_cons == NULL) {
      bgi->next_cons = bgi->cons_list;
    }

    /* now advance */
    while(bgi->next_cons != NULL) {
      cblk = bgi->next_cons->data;
      if(cblk->end >= p->pos) {
	      break;
      }
      bgi->next_cons = g_list_next(bgi->next_cons);
    }

    /* calculate B value and 1st/2nd derivatives at new position */
    bkgd_calc_b(p, bgi->cons_list, bgi->next_cons, bgi->parm, bgi->rtab->chr_r_len);

    /* calculate differences in position and b between points */
    b_delta = fabs(p->b - last_p->b);
    pos_delta = p->pos - last_p->pos;
    r_delta = (p->r_pos - last_p->r_pos);
   
    /*     fprintf(stderr, "B EVAL:p->pos=%ld, p->r_pos=%g,  p->b=%g, " */
    /* 	    "last_p->b=%g, pos_delta=%ld,r_delta=%g, b_delta=%g\n", */
    /*  	    p->pos, p->r_pos, p->b, last_p->b, pos_delta, r_delta, b_delta); */

    if (pos_delta > 1 && b_delta > BKGD_CHANGE_THRESH && r_delta > 0.0) {
      /* We don't want use this point because it is still too far
       * away.  Push it onto queue for later use.
       */
      g_queue_push_head(bgi->p_queue, bkgd_point_dup(p));
      
      /* Try a position that is 1/4 as far away as last position we
       * aimed for. Do not use r_delta here (which is always less than
       * or equal to the difference we were aiming for) because if
       * there is a jump in rec-dist at a particular site, we can end
       * up trying the same position over and over again.
       */
      // next_r = last_p->r_pos + (next_r - last_p->r_pos)*0.25;
      // DEBUG: try taking less conservative step 1/2 r_delta
      next_r = last_p->r_pos + (next_r - last_p->r_pos)*0.5;
    } 
    else {
      /* keep this position */
      keep_point = TRUE;
    }
  }

  /* update slope */
  if(r_delta == 0.0) {
    /* This occurs when two points are immediately adjacent (at least
     * in rec dist). In this case there is no change in the B values
     * between the points, so just call slope 0.
     */    
    bgi->slope = 0.0;
  } 
  else {
    bgi->slope = (p->b - last_p->b) / r_delta;
  }

}



/**
 * Creates a new BkgdInterp structure, that is intended to
 * transparently decide on appropriate locations to evaluate the B
 * function and to perform interpolation between the locations where B
 * was actually evaluated.
 */
// BkgdInterp *bkgd_interp_new(RecRateTable *rtab, long chr_len, 
//           GList *cons_list, BkgdParam *parm)
// DAVID: modify bgi init function to take Chromosome class
BkgdInterp *bkgd_interp_new(RecRateTable *rtab, Chromosome *chr, 
			    GList *cons_list, BkgdParam *parm) {
  BkgdInterp *bgi;

  bgi = g_new(BkgdInterp, 1);

  bgi->rtab = rtab;

  // DAVID: load Chromosome, set len from stored len
  bgi->chr = chr;
  bgi->chr_len = chr->len;
  bgi->cons_list = cons_list;
  bgi->next_cons = cons_list;
  assert (bgi->next_cons != NULL);

  // DAVID: last_cblk_end needs to be aligned to the current position using next cons
  bgi->last_cblk_rend = 0.0;

  // DAVID: align next_cons for each p_start
  ConsBlock *cblk = NULL;
  while (TRUE) {
    // no mores cons blocks: break
    if (bgi->next_cons == NULL)
      break;
    // next_cons points at a cblk: get data
    cblk = bgi->next_cons->data;
    // current block aligned: break
    if (cblk->end >= bgi->chr->p_start)
      break;
    // if this statement is reached next_cons must advance again: record last_cblk_rend first
    bgi->last_cblk_rend = rectab_rpos(bgi->rtab, cblk->end+1);
    // advance next_cons pointer
    bgi->next_cons = g_list_next(bgi->next_cons);
  }

  bgi->parm = parm;
  bgi->p_queue = g_queue_new();

  /* use first base on chr as first point */
  bgi->p1 = g_new(BkgdPoint, 1);
  // bgi->p1->pos = 1;
  // bgi->p1->r_pos = rectab_rpos(rtab, 1);
  // DAVID: use first base in PARTITION as first point
  bgi->p1->pos = bgi->chr->p_start;
  bgi->p1->r_pos = rectab_rpos(rtab, bgi->p1->pos);

  // calculate b, d1b, d2b at initial pos
  bkgd_calc_b(bgi->p1, bgi->cons_list, bgi->next_cons, parm, rtab->chr_r_len);

  bgi->p2 = g_new(BkgdPoint, 1);

  /* compute position of second point */
  get_next_point(bgi, bgi->p1, bgi->p2);

  return bgi;
}



/**
 * Frees memory allocated for BkgdInterp. 
 */
void bkgd_interp_free(BkgdInterp *bgi) {
  BkgdPoint *p;
  g_free(bgi->p1);
  g_free(bgi->p2);

  while(bgi->p_queue->length > 0) {
    p = g_queue_pop_head(bgi->p_queue);
    g_free(p);
  }
  g_queue_free(bgi->p_queue);

  g_free(bgi);
}


/**
 * Calculates an interpolated value at the desired position using the
 * provided initialised interpolator.
 */
// DEBUG: add param for debug file handle to write interpolator points, bvals to file
// double bkgd_interp_eval(BkgdInterp *bgi, long pos, FILE *d_fh) {
double bkgd_interp_eval(BkgdInterp *bgi, long pos) {

  double r_delta;
  BkgdPoint *p;

  if(bgi->p2->pos < pos) {
    /* need to find next point and update slope */
    p = bgi->p1;
    bgi->p1 = bgi->p2; /* second point becomes first point */
    get_next_point(bgi, bgi->p1, p);
    bgi->p2 = p; /* new point becomes second point */

    // DEBUG: write upper edge of each interpolator in bp and its exact b value
    // fprintf(d_fh, "%ld %.16f\n", bgi->p2->pos, bgi->p2->b);

    // DEBUG: write upper edge of each interpolator in bp and M and its exact b value
    // fprintf(d_fh, "%ld %.16g %.16g\n", bgi->p2->pos, bgi->p2->r_pos, bgi->p2->b);
    // fflush(d_fh);
    // DEBUG: write upper edge of interpolator in bp & M, exact b and bgi return condition
    // fprintf(d_fh, "%ld %.16g %.16g %s\n", bgi->p2->pos, bgi->p2->r_pos, bgi->p2->b, bgi->return_condition);
    // fflush(d_fh);

  }

  /* perform linear interpolation (in recombination space) */
  r_delta = rectab_rpos(bgi->rtab, pos) - bgi->p1->r_pos;

  /*   fprintf(stderr, "prev=%g, next=%g, slope=%g, r_delta=%g\n", */
  /* 	  bgi->p1->b, bgi->p2->b, bgi->slope, r_delta); */

  return bgi->p1->b + r_delta*bgi->slope;
}



