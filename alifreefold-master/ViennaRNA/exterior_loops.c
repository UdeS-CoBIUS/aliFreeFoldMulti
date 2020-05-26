#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/exterior_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

struct default_data {
  int                       *idx;
  unsigned char             *mx;
  unsigned char             **mx_window;
  int                       cp;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_callback_hc_evaluate *hc_f;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_el_t  *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ext_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      vrna_mx_pf_aux_el_t   *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ext_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           vrna_mx_pf_aux_el_t  *aux_mx);


PRIVATE int
E_ext_loop_5(vrna_fold_compound_t *vc);


PRIVATE int
E_ext_loop_5_comparative(vrna_fold_compound_t *vc);


PRIVATE int
BT_ext_loop_f5(vrna_fold_compound_t *fc,
               int                  *k,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count);


PRIVATE int
BT_ext_loop_f5_comparative(vrna_fold_compound_t *fc,
                           int                  *k,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count);


PRIVATE int
E_ext_loop_3(vrna_fold_compound_t *fc,
             int                  i);


PRIVATE int
E_ext_loop_3_comparative(vrna_fold_compound_t *fc,
                         int                  i);


PRIVATE int
BT_ext_loop_f3(vrna_fold_compound_t *vc,
               int                  *k,
               int                  maxdist,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count);


PRIVATE int
BT_ext_loop_f3_comparative(vrna_fold_compound_t *vc,
                           int                  *k,
                           int                  maxdist,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count);


PRIVATE int
BT_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  int                   *i,
                  int                   maxj);


PRIVATE int
BT_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              int                   *i,
                              int                   maxj);


PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data);


PRIVATE unsigned char
hc_default_window(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data);


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data);


PRIVATE unsigned char
hc_default_user_window(int            i,
                       int            j,
                       int            k,
                       int            l,
                       unsigned char  d,
                       void           *data);


PRIVATE INLINE int get_pair_type_md(int       i,
                                    int       j,
                                    vrna_md_t *md);


PRIVATE INLINE int
get_pair_type(int   ij,
              char  *ptype);


PRIVATE INLINE int
get_pair_type_window(int  i,
                     int  j,
                     char **ptype);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_ext_loop_5(vrna_fold_compound_t *fc)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return E_ext_loop_5(fc);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return E_ext_loop_5_comparative(fc);
        break;
    }
  }

  return 0;
}


PUBLIC int
vrna_BT_ext_loop_f5(vrna_fold_compound_t  *fc,
                    int                   *k,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    int                   *stack_count)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return BT_ext_loop_f5(fc, k, i, j, bp_stack, stack_count);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return BT_ext_loop_f5_comparative(fc, k, i, j, bp_stack, stack_count);
        break;
    }
  }

  return -1;
}


PUBLIC int
vrna_E_ext_loop_3(vrna_fold_compound_t  *fc,
                  int                   i)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return E_ext_loop_3(fc, i);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return E_ext_loop_3_comparative(fc, i);
        break;
    }
  }

  return 0;
}


PUBLIC int
vrna_BT_ext_loop_f3(vrna_fold_compound_t  *fc,
                    int                   *k,
                    int                   maxdist,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    int                   *stack_count)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return BT_ext_loop_f3(fc, k, maxdist, i, j, bp_stack, stack_count);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return BT_ext_loop_f3_comparative(fc, k, maxdist, i, j, bp_stack, stack_count);
        break;
    }
  }

  return -1;
}


PUBLIC int
vrna_BT_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       int                  *i,
                       int                  maxj)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return BT_ext_loop_f3_pp(fc, i, maxj);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return BT_ext_loop_f3_pp_comparative(fc, i, maxj);
        break;
    }
  }

  return -1;
}


PUBLIC int
E_Stem(int          type,
       int          si1,
       int          sj1,
       int          extLoop,
       vrna_param_t *P)
{
  int energy  = 0;
  int d5      = (si1 >= 0) ? P->dangle5[type][si1] : 0;
  int d3      = (sj1 >= 0) ? P->dangle3[type][sj1] : 0;

  if (type > 2)
    energy += P->TerminalAU;

  if (si1 >= 0 && sj1 >= 0)
    energy += (extLoop) ? P->mismatchExt[type][si1][sj1] : P->mismatchM[type][si1][sj1];
  else
    energy += d5 + d3;

  if (!extLoop)
    energy += P->MLintern[type];

  return energy;
}


PUBLIC int
E_ExtLoop(int           type,
          int           si1,
          int           sj1,
          vrna_param_t  *P)
{
  int energy = 0;

  if (si1 >= 0 && sj1 >= 0)
    energy += P->mismatchExt[type][si1][sj1];
  else if (si1 >= 0)
    energy += P->dangle5[type][si1];
  else if (sj1 >= 0)
    energy += P->dangle3[type][sj1];

  if (type > 2)
    energy += P->TerminalAU;

  return energy;
}


PUBLIC FLT_OR_DBL
exp_E_Stem(int              type,
           int              si1,
           int              sj1,
           int              extLoop,
           vrna_exp_param_t *P)
{
  double  energy  = 1.0;
  double  d5      = (si1 >= 0) ? P->expdangle5[type][si1] : 1.;
  double  d3      = (sj1 >= 0) ? P->expdangle3[type][sj1] : 1.;

  if (si1 >= 0 && sj1 >= 0)
    energy = (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy = d5 * d3;

  if (type > 2)
    energy *= P->expTermAU;

  if (!extLoop)
    energy *= P->expMLintern[type];

  return (FLT_OR_DBL)energy;
}


PUBLIC FLT_OR_DBL
exp_E_ExtLoop(int               type,
              int               si1,
              int               sj1,
              vrna_exp_param_t  *P)
{
  double energy = 1.0;

  if (si1 >= 0 && sj1 >= 0)
    energy = P->expmismatchExt[type][si1][sj1];
  else if (si1 >= 0)
    energy = P->expdangle5[type][si1];
  else if (sj1 >= 0)
    energy = P->expdangle3[type][sj1];

  if (type > 2)
    energy *= P->expTermAU;

  return (FLT_OR_DBL)energy;
}


PUBLIC int
E_ext_loop(int                  i,
           int                  j,
           vrna_fold_compound_t *vc)
{
  char                      *ptype;
  unsigned char             *hard_constraints;
  short                     *S;
  int                       ij, en, e, type, cp, *idx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp                = vc->cutpoint;
  S                 = vc->sequence_encoding;
  idx               = vc->jindx;
  ptype             = vc->ptype;
  P                 = vc->params;
  md                = &(P->model_details);
  hard_constraints  = vc->hc->matrix;
  sc                = vc->sc;

  hc_dat_local.idx  = idx;
  hc_dat_local.mx   = hard_constraints;
  hc_dat_local.cp   = cp;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  e     = INF;
  ij    = idx[j] + i;
  type  = get_pair_type(ij, ptype);

  if ((cp < 0) || (((i) >= cp) || ((j) < cp))) {
    /* regular exterior loop */
    if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
      switch (md->dangles) {
        case 2:
          e = E_ExtLoop(type, S[i - 1], S[j + 1], P);
          break;

        case 0:
        /* fall through */

        default:
          e = E_ExtLoop(type, -1, -1, P);
          break;
      }
      if (sc)
        if (sc->f)
          e += sc->f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
    }

    if (md->dangles % 2) {
      ij = idx[j - 1] + i;
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = get_pair_type(ij, ptype);

        en = E_ExtLoop(type, -1, S[j], P);

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);

        e = MIN2(e, en);
      }

      ij = idx[j] + i + 1;
      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = get_pair_type(ij, ptype);

        en = E_ExtLoop(type, S[i], -1, P);

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

        e = MIN2(e, en);
      }
    }
  }

  return e;
}


PRIVATE INLINE int
get_pair_type_md(int        i,
                 int        j,
                 vrna_md_t  *md)
{
  int tt = md->pair[i][j];

  return (tt == 0) ? 7 : tt;
}


PRIVATE INLINE int
get_pair_type(int   ij,
              char  *ptype)
{
  int tt = (int)ptype[ij];

  return (tt == 0) ? 7 : tt;
}


PRIVATE INLINE int
get_pair_type_window(int  i,
                     int  j,
                     char **ptype)
{
  int tt = (int)ptype[i][j - i];

  return (tt == 0) ? 7 : tt;
}


PRIVATE int
E_ext_loop_5(vrna_fold_compound_t *vc)
{
  char                      *ptype;
  unsigned char             *hc;
  short                     *S;
  int                       en, i, j, ij, type, length, *indx, *hc_up, *f5, *c, dangle_model,
                            *ggg, with_gquad, turn, k, u, with_ud;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = (int)vc->length;
  ptype         = vc->ptype;
  S             = vc->sequence_encoding;
  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_ext;
  sc            = vc->sc;
  f5            = vc->matrices->f5;
  c             = vc->matrices->c;
  P             = vc->params;
  dangle_model  = P->model_details.dangles;
  ggg           = vc->matrices->ggg;
  with_gquad    = P->model_details.gquad;
  turn          = P->model_details.min_loop_size;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;

  hc_dat_local.idx    = indx;
  hc_dat_local.mx     = hc;
  hc_dat_local.hc_up  = hc_up;
  hc_dat_local.cp     = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  f5[0] = 0;
  for (i = 1; i <= turn + 1; i++) {
    if (f5[i - 1] != INF) {
      if (evaluate(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        f5[i] = f5[i - 1];

        if (sc) {
          if (sc->energy_up)
            f5[i] += sc->energy_up[i][1];

          if (sc->f)
            f5[i] += sc->f(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }
      } else {
        f5[i] = INF;
      }
    } else {
      f5[i] = INF;
    }
  }

  if (with_ud) {
    /* do we include ligand binding? */
    /*  construct all possible combinations of
     *  f[i-1] + L[i,j] with j <= turn + 1
     */
    for (i = 1; i <= turn + 1; i++) {
      if (f5[i - 1] != INF) {
        for (k = 0; k < domains_up->uniq_motif_count; k++) {
          u = domains_up->uniq_motif_size[k];
          j = i + u - 1;
          if (j <= turn + 1) {
            if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
              en = f5[i - 1] +
                   domains_up->energy_cb(vc,
                                         i, j,
                                         VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                         domains_up->data);

              if (sc) {
                if (sc->energy_up)
                  en += sc->energy_up[i][u];

                if (sc->f)
                  en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
              }

              f5[j] = MIN2(f5[j], en);
            }
          }
        }
      }
    }
  }

  /* duplicated code may be faster than conditions inside loop ;) */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = turn + 2; j <= length; j++) {
        /* initialize with INF */
        f5[j] = INF;

        /* check for 3' extension with one unpaired nucleotide */
        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (sc) {
              if (sc->energy_up)
                f5[j] += sc->energy_up[j][1];

              if (sc->f)
                f5[j] += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
            }
          }
        }

        if (with_ud) {
          for (k = 0; k < domains_up->uniq_motif_count; k++) {
            u = domains_up->uniq_motif_size[k];
            if ((j - u >= 0) && (f5[j - u] != INF)) {
              if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
                en = f5[j - u] +
                     domains_up->energy_cb(vc,
                                           j - u + 1,
                                           j,
                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                           domains_up->data);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j - u + 1][u];

                  if (sc->f)
                    en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        /* check for possible stems branching off the exterior loop */
        if (sc && sc->f) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, -1, -1, P) +
                       sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, -1, -1, P);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, -1, P);

            if (sc)
              if (sc->f)
                en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

            f5[j] = MIN2(f5[j], en);
          }
        }
      }
      break;

    /* always use dangles on both sides */
    case 2:
      for (j = turn + 2; j < length; j++) {
        f5[j] = INF;

        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (sc) {
              if (sc->energy_up)
                f5[j] += sc->energy_up[j][1];

              if (sc->f)
                f5[j] += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
            }
          }
        }

        if (with_ud) {
          for (k = 0; k < domains_up->uniq_motif_count; k++) {
            u = domains_up->uniq_motif_size[k];
            if ((j - u >= 0) && (f5[j - u] != INF)) {
              if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
                en = f5[j - u] +
                     domains_up->energy_cb(vc,
                                           j - u + 1,
                                           j,
                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                           domains_up->data);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j - u + 1][u];

                  if (sc->f)
                    en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        if (sc && sc->f) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, S[i - 1], S[j + 1], P) +
                       sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, S[i - 1], S[j + 1], P);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, S[j + 1], P);

            if (sc)
              if (sc->f)
                en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

            f5[j] = MIN2(f5[j], en);
          }
        }
      }

      f5[length] = INF;
      if (f5[length - 1] != INF) {
        if (evaluate(1, length, 1, length - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          f5[length] = f5[length - 1];
          if (sc) {
            if (sc->energy_up)
              f5[length] += sc->energy_up[length][1];

            if (sc->f)
              f5[length] += sc->f(1, length, 1, length - 1, VRNA_DECOMP_EXT_EXT, sc->data);
          }
        }
      }

      if (with_ud) {
        for (k = 0; k < domains_up->uniq_motif_count; k++) {
          u = domains_up->uniq_motif_size[k];
          if ((length - u >= 0) && (f5[length - u] != INF)) {
            if (evaluate(1, length, 1, length - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
              en = f5[length - u] +
                   domains_up->energy_cb(vc,
                                         length - u + 1,
                                         length,
                                         VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);

              if (sc) {
                if (sc->energy_up)
                  en += sc->energy_up[length - u + 1][u];

                if (sc->f)
                  en += sc->f(1, length, 1, length - u, VRNA_DECOMP_EXT_EXT, sc->data);
              }

              f5[length] = MIN2(f5[length], en);
            }
          }
        }
      }

      if (sc && sc->f) {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, S[i - 1], -1, P) +
                     sc->f(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      } else {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, S[i - 1], -1, P);

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      }

      ij = indx[length] + 1;

      if (with_gquad)
        f5[length] = MIN2(f5[length], ggg[ij]);

      if (c[ij] != INF) {
        if (evaluate(1, length, 1, length, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          type = get_pair_type(ij, ptype);

          en = c[ij] +
               E_ExtLoop(type, -1, -1, P);

          if (sc)
            if (sc->f)
              en += sc->f(1, length, 1, length, VRNA_DECOMP_EXT_STEM, sc->data);

          f5[length] = MIN2(f5[length], en);
        }
      }

      break;

    /* normal dangles, aka dangle_model = 1 || 3 */
    default:
      for (j = turn + 2; j <= length; j++) {
        f5[j] = INF;
        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (sc) {
              if (sc->energy_up)
                f5[j] += sc->energy_up[j][1];

              if (sc->f)
                f5[j] += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
            }
          }
        }

        if (with_ud) {
          for (k = 0; k < domains_up->uniq_motif_count; k++) {
            u = domains_up->uniq_motif_size[k];
            if ((j - u >= 0) && (f5[j - u] != INF)) {
              if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
                en = f5[j - u] +
                     domains_up->energy_cb(vc,
                                           j - u + 1,
                                           j,
                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                           domains_up->data);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j - u + 1][u];

                  if (sc->f)
                    en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        for (i = j - turn - 1; i > 1; i--) {
          ij = indx[j] + i;
          if (f5[i - 1] != INF) {
            if (with_gquad)
              f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, -1, -1, P);

                if (sc)
                  if (sc->f)
                    en += sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                f5[j] = MIN2(f5[j], en);
              }
            }
          }

          if ((f5[i - 2] != INF) && c[ij] != INF) {
            if (evaluate(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
              type = get_pair_type(ij, ptype);

              en = f5[i - 2] + c[ij] +
                   E_ExtLoop(type, S[i - 1], -1, P);

              if (sc) {
                if (sc->energy_up)
                  en += sc->energy_up[i - 1][1];

                if (sc->f)
                  en += sc->f(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
              }

              f5[j] = MIN2(f5[j], en);
            }
          }

          ij = indx[j - 1] + i;
          if (c[ij] != INF) {
            if (f5[i - 1] != INF) {
              if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, -1, S[j], P);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j][1];

                  if (sc->f)
                    en += sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }

            if (f5[i - 2] != INF) {
              if (evaluate(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 2] +
                     c[ij] +
                     E_ExtLoop(type, S[i - 1], S[j], P);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[i - 1][1] +
                          sc->energy_up[j][1];

                  if (sc->f)
                    en += sc->f(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, -1, P);

            if (sc)
              if (sc->f)
                en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

            f5[j] = MIN2(f5[j], en);
          }
        }

        ij = indx[j - 1] + 1;
        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, S[j], P);

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[j][1];

              if (sc->f)
                en += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            f5[j] = MIN2(f5[j], en);
          }
        }
      }         /* end for j... */
      break;
  }

  return f5[length];
}


PRIVATE int
E_ext_loop_5_comparative(vrna_fold_compound_t *vc)
{
  unsigned char             *hc;
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       en, i, j, ij, tt, length, *indx, *hc_up, *f5, *c, dangle_model,
                            *ggg, with_gquad, turn, n_seq, s;
  vrna_sc_t                 **scs;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = (int)vc->length;
  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;
  S3            = vc->S3;
  a2s           = vc->a2s;
  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_ext;
  scs           = vc->scs;
  f5            = vc->matrices->f5;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;

  hc_dat_local.idx    = indx;
  hc_dat_local.mx     = hc;
  hc_dat_local.hc_up  = hc_up;
  hc_dat_local.cp     = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  f5[0] = 0;
  for (i = 1; i <= turn + 1; i++) {
    if (f5[i - 1] != INF) {
      if (evaluate(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        f5[i] = f5[i - 1];

        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->energy_up)
                f5[i] += scs[s]->energy_up[a2s[s][j]][1];

              if (scs[s]->f)
                f5[i] += scs[s]->f(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
            }
          }
        }
      } else {
        f5[i] = INF;
      }
    } else {
      f5[i] = INF;
    }
  }

  /* duplicated code may be faster than conditions inside loop ;) */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = turn + 2; j <= length; j++) {
        /* initialize with INF */
        f5[j] = INF;

        /* check for 3' extension with one unpaired nucleotide */
        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->energy_up)
                    f5[j] += scs[s]->energy_up[a2s[s][j]][1];

                  if (scs[s]->f)
                    f5[j] += scs[s]->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
                }
              }
            }
          }
        }

        /* check for possible stems branching off the exterior loop */
        if (scs) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] + c[ij];
                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    en  += E_ExtLoop(tt, -1, -1, P);

                    if (scs[s] && scs[s]->f)
                      en += scs[s]->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, scs[s]->data);
                  }
                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] +
                       c[ij];

                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    en  += E_ExtLoop(tt, -1, -1, P);
                  }

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            en = c[ij];
            if (scs) {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                en  += E_ExtLoop(tt, -1, -1, P);

                if (scs[s] && scs[s]->f)
                  en += scs[s]->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
              }
            } else {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                en  += E_ExtLoop(tt, -1, -1, P);
              }
            }

            f5[j] = MIN2(f5[j], en);
          }
        }
      }
      break;

    /* always use dangles on both sides */
    case 2:
      for (j = turn + 2; j < length; j++) {
        f5[j] = INF;

        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->energy_up)
                    f5[j] += scs[s]->energy_up[a2s[s][j]][1];

                  if (scs[s]->f)
                    f5[j] += scs[s]->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
                }
              }
            }
          }
        }

        if (scs) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] +
                       c[ij];

                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    en  += E_ExtLoop(tt, S5[s][i], S3[s][j], P);

                    if (scs[s] && scs[s]->f)
                      en += scs[s]->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, scs[s]->data);
                  }
                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] +
                       c[ij];

                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    en  += E_ExtLoop(tt, S5[s][i], S3[s][j], P);
                  }
                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            en = c[ij];

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                en  += E_ExtLoop(tt, -1, S3[s][j], P);

                if (scs[s] && scs[s]->f)
                  en += scs[s]->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
              }
            } else {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                en  += E_ExtLoop(tt, -1, S3[s][j], P);
              }
            }

            f5[j] = MIN2(f5[j], en);
          }
        }
      }

      f5[length] = INF;
      if (f5[length - 1] != INF) {
        if (evaluate(1, length, 1, length - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          f5[length] = f5[length - 1];
          if (scs) {
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                if (scs[s]->energy_up)
                  f5[length] += scs[s]->energy_up[a2s[s][length]][1];

                if (scs[s]->f) {
                  f5[length] += scs[s]->f(1,
                                          length,
                                          1,
                                          length - 1,
                                          VRNA_DECOMP_EXT_EXT,
                                          scs[s]->data);
                }
              }
            }
          }
        }
      }

      if (scs) {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                en = f5[i - 1] +
                     c[ij];

                for (s = 0; s < n_seq; s++) {
                  tt  = get_pair_type_md(S[s][i], S[s][j], md);
                  en  += E_ExtLoop(tt, S5[s][i], -1, P);

                  if (scs[s] && scs[s]->f)
                    en += scs[s]->f(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, scs[s]->data);
                }

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      } else {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                en = f5[i - 1] +
                     c[ij];

                for (s = 0; s < n_seq; s++) {
                  tt  = get_pair_type_md(S[s][i], S[s][j], md);
                  en  += E_ExtLoop(tt, S5[s][i], -1, P);
                }

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      }

      ij = indx[length] + 1;

      if (with_gquad)
        f5[length] = MIN2(f5[length], ggg[ij]);

      if (c[ij] != INF) {
        if (evaluate(1, length, 1, length, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          en = c[ij];

          if (scs) {
            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][i], S[s][j], md);
              en  += E_ExtLoop(tt, -1, -1, P);

              if (scs[s] && scs[s]->f)
                en += scs[s]->f(1, length, 1, length, VRNA_DECOMP_EXT_STEM, scs[s]->data);
            }
          } else {
            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][i], S[s][j], md);
              en  += E_ExtLoop(tt, -1, -1, P);
            }
          }

          f5[length] = MIN2(f5[length], en);
        }
      }

      break;
  }

  return f5[length];
}


PRIVATE int
E_ext_loop_3(vrna_fold_compound_t *fc,
             int                  i)
{
  char                      **ptype;
  short                     *S1;
  int                       e, dangle_model, *f3, j, turn, length, maxdist, with_gquad, **ggg,
                            energy, type, **c;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  e = INF;

  length        = fc->length;
  maxdist       = fc->window_size;
  S1            = fc->sequence_encoding;
  ptype         = fc->ptype_local;
  P             = fc->params;
  md            = &(P->model_details);
  hc            = fc->hc;
  sc            = fc->sc;
  f3            = fc->matrices->f3_local;
  c             = fc->matrices->c_local;
  ggg           = fc->matrices->ggg_local;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;

  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.cp         = fc->cutpoint;

  if (fc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = fc->hc->f;
    hc_dat_local.hc_dat = fc->hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /* first case: i stays unpaired */
  if (evaluate(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    e = f3[i + 1];
    if (sc) {
      if (sc->energy_up)
        e += sc->energy_up[i][1];

      if (sc->f)
        e += sc->f(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);
    }
  }

  /* next all cases where i is paired */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (f3[j + 1] != INF) && (ggg[i][j - i] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            type = get_pair_type_window(i, j, ptype);

            energy = f3[j + 1] +
                     c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, length, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            type = get_pair_type_window(i, j, ptype);

            energy = c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

            e = MIN2(e, energy);
          }
        }
      }

      break;
    /* always use dangle_model on both sides */
    case 2:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (ggg[i][j - i] != INF) && (f3[j + 1] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            type = get_pair_type_window(i, j, ptype);

            energy = f3[j + 1] +
                     c[i][j - i] +
                     E_ExtLoop(type,
                               (i > 1) ? S1[i - 1] : -1,
                               S1[j + 1],
                               P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            type = get_pair_type_window(i, j, ptype);

            energy = c[i][j - i] +
                     E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

            e = MIN2(e, energy);
          }
        }
      }

      break;
    /* normal dangle_model, aka dangle_model = 1 */
    default:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if (with_gquad && (f3[j + 1] != INF) && (ggg[i][j - i] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        type = get_pair_type_window(i, j, ptype);

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            energy = f3[j + 1] +
                     c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            e = MIN2(e, energy);
          }
        }

        if (j + 2 <= length) {
          if (evaluate(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            if ((c[i][j - i] != INF) && (f3[j + 2] != INF)) {
              energy = c[i][j - i] +
                       f3[j + 2] +
                       E_ExtLoop(type, -1, S1[j + 1], P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[j + 1][1];

                if (sc->f)
                  energy += sc->f(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
              }

              e = MIN2(e, energy);
            }
          }
        } else {
          if (evaluate(i, length, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            if (c[i][j - i] != INF) {
              energy = c[i][j - i] +
                       E_ExtLoop(type, -1, S1[j + 1], P);

              if ((sc) && (sc->f))
                energy += sc->f(i, length, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

              e = MIN2(e, energy);
            }
          }
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
          type = get_pair_type_window(i + 1, j, ptype);

          if ((c[i + 1][j - i - 1] != INF) && (f3[j + 1] != INF)) {
            energy = f3[j + 1] +
                     c[i + 1][j - i - 1] +
                     E_ExtLoop(type, S1[i], -1, P);

            if (sc) {
              if (sc->energy_up)
                energy += sc->energy_up[i][1];

              if (sc->f)
                energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
            }

            e = MIN2(e, energy);
          }
        }

        if (j + 2 <= length) {
          if (evaluate(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            if ((c[i + 1][j - i - 1] != INF) && (f3[j + 2] != INF)) {
              energy = c[i + 1][j - i - 1] +
                       f3[j + 2] +
                       E_ExtLoop(type, S1[i], S1[j + 1], P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i][1] +
                            sc->energy_up[j + 1][1];

                if (sc->f)
                  energy += sc->f(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
              }

              e = MIN2(e, energy);
            }
          }
        } else {
          if (evaluate(i, length, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            if (c[i + 1][j - i - 1] != INF) {
              energy = c[i + 1][j - i - 1] +
                       E_ExtLoop(type, S1[i], S1[j + 1], P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i][1];

                if (sc->f)
                  energy += sc->f(i, length, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);
              }

              e = MIN2(e, energy);
            }
          }
        }
      }

      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            type = get_pair_type_window(i, j, ptype);

            energy = c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

            e = MIN2(e, energy);
          }
        }

        if (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i + 1][j - i - 1] != INF) {
            type = get_pair_type_window(i + 1, j, ptype);

            energy = c[i + 1][j - i - 1] +
                     E_ExtLoop(type, S1[i], -1, P);

            if (sc) {
              if (sc->energy_up)
                energy += sc->energy_up[i][1];

              if (sc->f)
                energy += sc->f(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            e = MIN2(e, energy);
          }
        }
      }

      break;
  } /* switch(dangle_model)... */

  return e;
}


PRIVATE int
E_ext_loop_3_comparative(vrna_fold_compound_t *fc,
                         int                  i)
{
  short                     **S, **S5, **S3;
  int                       e, dangle_model, *f3, j, turn, length, maxdist, with_gquad, **ggg,
                            energy, **c, n_seq, s, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  e = INF;

  length        = fc->length;
  n_seq         = fc->n_seq;
  S             = fc->S;
  S5            = fc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
  S3            = fc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
  maxdist       = fc->window_size;
  P             = fc->params;
  md            = &(P->model_details);
  hc            = fc->hc;
  f3            = fc->matrices->f3_local;
  c             = fc->matrices->c_local;
  ggg           = fc->matrices->ggg_local;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;

  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.cp         = fc->cutpoint;

  if (fc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = fc->hc->f;
    hc_dat_local.hc_dat = fc->hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /* first case: i stays unpaired */
  if (evaluate(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local))
    e = f3[i + 1];

  /* next all cases where i is paired */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (f3[j + 1] != INF) && (ggg[i][j - i] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            energy = f3[j + 1] +
                     c[i][j - i];

            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, -1, -1, P);
            }

            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, length, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            energy = c[i][j - i];
            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, -1, -1, P);
            }
            e = MIN2(e, energy);
          }
        }
      }

      break;
    /* always use dangle_model on both sides */
    case 2:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (ggg[i][j - i] != INF) && (f3[j + 1] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            energy = f3[j + 1] +
                     c[i][j - i];

            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, (i > 1) ? S5[s][i] : -1, S3[s][j], P);
            }
            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            energy = c[i][j - i];
            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, (i > 1) ? S5[s][i] : -1, -1, P);
            }
            e = MIN2(e, energy);
          }
        }
      }

      break;
  } /* switch(dangle_model)... */

  return e;
}


PRIVATE int
BT_ext_loop_f5(vrna_fold_compound_t *vc,
               int                  *k,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count)
{
  char                      *ptype;
  short                     mm5, mm3, *S1;
  unsigned int              *sn;
  int                       length, fij, fi, jj, u, en, e, *my_f5, *my_c, *my_ggg, *idx,
                            dangle_model, turn, with_gquad, cnt, ii, with_ud, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = vc->length;
  P             = vc->params;
  md            = &(P->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  my_f5         = vc->matrices->f5;
  my_c          = vc->matrices->c;
  my_ggg        = vc->matrices->ggg;
  domains_up    = vc->domains_up;
  idx           = vc->jindx;
  ptype         = vc->ptype;
  S1            = vc->sequence_encoding;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;

  hc_dat_local.idx    = idx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  jj = *k;

  /* nibble off unpaired 3' stretches harboring bound ligands (interspersed with unpaired nucleotides) */
  if (with_ud) {
    do {
      fij = my_f5[jj];
      fi  = INF;

      /* try nibble off one unpaired nucleotide first */
      if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_f5[jj - 1];

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (jj == 1) {
          /* no more pairs */
          *i  = *j = -1;
          *k  = 0;
          return 1;
        }

        if (fij == fi) {
          jj--;
          continue;
        }
      }

      /* next, try nibble off a ligand */
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        ii  = jj - u + 1;
        if ((ii > 0) && evaluate(1, jj, 1, jj - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          en = domains_up->energy_cb(vc,
                                     ii,
                                     jj,
                                     VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[ii][u];

            if (sc->f)
              en += sc->f(1, jj, 1, jj - u, VRNA_DECOMP_EXT_EXT, sc->data);
          }

          fi  = my_f5[ii - 1];
          fi  += en;

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            jj = ii - 1;
            break;
          }
        }
      }

      if (jj == 0) {
        /* no more pairs */
        *i  = *j = -1;
        *k  = 0;
        return 1;
      }
    } while (fij == fi);
  } else {
    /* nibble off unpaired 3' bases */
    do {
      fij = my_f5[jj];
      fi  = INF;

      if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_f5[jj - 1];

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }
      }

      if (--jj == 0)
        break;
    } while (fij == fi);
    jj++;
  }

  if (jj < turn + 2) {
    /* no more pairs */
    *i  = *j = -1;
    *k  = 0;
    return 1;
  }

  /* must have found a decomposition */
  switch (dangle_model) {
    case 0:   /* j is paired. Find pairing partner */
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          type = get_pair_type(idx[jj] + u, ptype);

          en = my_c[idx[jj] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          if (sn[jj] != sn[u])
            en += P->DuplexInit;

          if (fij == E_ExtLoop(type, -1, -1, P) + en + my_f5[u - 1]) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;

    case 2:
      mm3 = ((jj < length) && (sn[jj + 1] == sn[jj])) ? S1[jj + 1] : -1;
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          mm5   = ((u > 1) && (sn[u] == sn[u - 1])) ? S1[u - 1] : -1;
          type  = get_pair_type(idx[jj] + u, ptype);

          en = my_c[idx[jj] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          if (sn[jj] != sn[u])
            en += P->DuplexInit;

          if (fij == E_ExtLoop(type, mm5, mm3, P) + en + my_f5[u - 1]) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;

    default:
      if (with_gquad) {
        if (fij == my_ggg[idx[jj] + 1]) {
          *i  = *j = -1;
          *k  = 0;
          vrna_BT_gquad_mfe(vc, 1, jj, bp_stack, stack_count);
          return 1;
        }
      }

      if (evaluate(1, jj, 1, jj, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = get_pair_type(idx[jj] + 1, ptype);

        en = my_c[idx[jj] + 1];
        if (sc)
          if (sc->f)
            en += sc->f(1, jj, 1, jj, VRNA_DECOMP_EXT_STEM, sc->data);

        if (sn[jj] != sn[1])
          en += P->DuplexInit;

        if (fij == en + E_ExtLoop(type, -1, -1, P)) {
          *i                            = 1;
          *j                            = jj;
          *k                            = 0;
          bp_stack[++(*stack_count)].i  = 1;
          bp_stack[(*stack_count)].j    = jj;
          return 1;
        }
      }

      if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        if (sn[jj] == sn[jj - 1]) {
          mm3   = S1[jj];
          type  = get_pair_type(idx[jj - 1] + 1, ptype);

          en = my_c[idx[jj - 1] + 1];
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[jj][1];

            if (sc->f)
              en += sc->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_STEM, sc->data);
          }

          if (sn[jj - 1] != sn[1])
            en += P->DuplexInit;

          if (fij == en + E_ExtLoop(type, -1, mm3, P)) {
            *i                            = 1;
            *j                            = jj - 1;
            *k                            = 0;
            bp_stack[++(*stack_count)].i  = 1;
            bp_stack[(*stack_count)].j    = jj - 1;
            return 1;
          }
        }
      }

      for (u = jj - turn - 1; u > 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        type = get_pair_type(idx[jj] + u, ptype);

        en = my_c[idx[jj] + u];
        if (sn[jj] != sn[u])
          en += P->DuplexInit;

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          e = my_f5[u - 1] +
              en +
              E_ExtLoop(type, -1, -1, P);

          if (sc)
            if (sc->f)
              e += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          if (fij == e) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }

        if (evaluate(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          if (sn[u] == sn[u - 1]) {
            mm5 = S1[u - 1];
            e   = my_f5[u - 2] +
                  en +
                  E_ExtLoop(type, mm5, -1, P);

            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[u - 1][1];

              if (sc->f)
                e += sc->f(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
            }

            if (fij == e) {
              *i                            = u;
              *j                            = jj;
              *k                            = u - 2;
              bp_stack[++(*stack_count)].i  = u;
              bp_stack[(*stack_count)].j    = jj;
              return 1;
            }
          }
        }

        type = get_pair_type(idx[jj - 1] + u, ptype);

        en = my_c[idx[jj - 1] + u];
        if (sn[jj - 1] != sn[u])
          en += P->DuplexInit;

        mm5 = (sn[u] == sn[u - 1]) ? S1[u - 1] : -1;
        mm3 = (sn[jj] == sn[jj - 1]) ? S1[jj] : -1;

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
          e = my_f5[u - 1] +
              en +
              E_ExtLoop(type, -1, mm3, P);

          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[jj][1];

            if (sc->f)
              e += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
          }

          if (fij == e) {
            *i                            = u;
            *j                            = jj - 1;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj - 1;
            return 1;
          }
        }

        if (evaluate(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
          e = my_f5[u - 2] + en + E_ExtLoop(type, mm5, mm3, P);
          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[jj][1] +
                   sc->energy_up[u - 1][1];

            if (sc->f)
              e += sc->f(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
          }

          if (fij == e) {
            *i                            = u;
            *j                            = jj - 1;
            *k                            = u - 2;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj - 1;
            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f5_comparative(vrna_fold_compound_t *vc,
                           int                  *k,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count)
{
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       length, fij, fi, jj, u, en, *my_f5, *my_c, *my_ggg, *idx,
                            dangle_model, turn, with_gquad, n_seq, ss, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = vc->length;
  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;
  S3            = vc->S3;
  a2s           = vc->a2s;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  my_f5         = vc->matrices->f5;
  my_c          = vc->matrices->c;
  my_ggg        = vc->matrices->ggg;
  idx           = vc->jindx;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;

  hc_dat_local.idx    = idx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  jj = *k;

  /* nibble off unpaired 3' bases */
  do {
    fij = my_f5[jj];
    fi  = INF;

    if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fi = my_f5[jj - 1];

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[a2s[ss][jj]][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, scs[ss]->data);
          }
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  if (jj < turn + 2) {
    /* no more pairs */
    *i  = *j = -1;
    *k  = 0;
    return 1;
  }

  /* must have found a decomposition */
  switch (dangle_model) {
    case 0:   /* j is paired. Find pairing partner */
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          en = my_c[idx[jj] + u] +
               my_f5[u - 1];

          for (ss = 0; ss < n_seq; ss++) {
            tt  = get_pair_type_md(S[ss][u], S[ss][jj], md);
            en  += E_ExtLoop(tt, -1, -1, P);
          }

          if (scs) {
            for (ss = 0; ss < n_seq; ss++)
              if (scs[ss] && scs[ss]->f)
                en += scs[ss]->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, scs[ss]->data);
          }

          if (fij == en) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;

    case 2:
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          en = my_c[idx[jj] + u] +
               my_f5[u - 1];

          for (ss = 0; ss < n_seq; ss++) {
            tt  = get_pair_type_md(S[ss][u], S[ss][jj], md);
            en  += E_ExtLoop(tt, (u > 1) ? S5[ss][u] : -1, (jj < length) ? S3[ss][jj] : -1, P);
          }

          if (scs) {
            for (ss = 0; ss < n_seq; ss++)
              if (scs[ss] && scs[ss]->f)
                en += scs[ss]->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, scs[ss]->data);
          }

          if (fij == en) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f3(vrna_fold_compound_t *vc,
               int                  *k,
               int                  maxdist,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count)
{
  char                      **ptype;
  short                     mm5, mm3, *S1;
  int                       length, fij, fj, ii, u, *f3, **c, **ggg, *idx,
                            dangle_model, turn, with_gquad, en, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = vc->length;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  sc            = vc->sc;
  f3            = vc->matrices->f3_local;
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  idx           = vc->jindx;
  ptype         = vc->ptype_local;
  S1            = vc->sequence_encoding;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;

  hc_dat_local.idx        = idx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.cp         = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  ii = *k;

  /* nibble off unpaired 5' bases */
  do {
    fij = f3[ii];
    fj  = INF;

    if (evaluate(ii, length, ii + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fj = f3[ii + 1];
      if (sc) {
        if (sc->energy_up)
          fj += sc->energy_up[ii][1];

        if (sc->f)
          fj += sc->f(ii, length, ii + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);
      }
    }

    if (++ii > maxdist)
      break;
  } while (fij == fj);
  ii--;

  if (ii > maxdist - turn + 1) {
    /* no more pairs */
    *i  = *j = -1;
    *k  = 0;
    return 1;
  }

  /*
   *  must have found a decomposition
   *  i is paired. Find pairing partner
   */
  switch (dangle_model) {
    /* no dangles */
    case 0:
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            vrna_BT_gquad_mfe(vc, ii, u, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          type = get_pair_type_window(ii, u, ptype);

          en = c[ii][u - ii] +
               E_ExtLoop(type, -1, -1, P) +
               f3[u + 1];

          if ((sc) && (sc->f))
            en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      mm5 = (ii > 1) ? S1[ii - 1] : -1;
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            vrna_BT_gquad_mfe(vc, ii, u, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          mm3   = (u < length) ? S1[u + 1] : -1;
          type  = get_pair_type_window(ii, u, ptype);

          en = c[ii][u - ii] +
               E_ExtLoop(type, mm5, mm3, P) +
               f3[u + 1];

          if ((sc) && (sc->f))
            en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    default:
      mm5 = S1[ii];
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            vrna_BT_gquad_mfe(vc, ii, u, bp_stack, stack_count);
            return 1;
          }
        }

        if (u + 2 <= length) {
          if (evaluate(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            mm3   = S1[u + 1];
            type  = get_pair_type_window(ii + 1, u, ptype);

            en = c[ii + 1][u - ii - 1] + E_ExtLoop(type, mm5, mm3, P) + f3[u + 2];

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[u + 1][1] +
                      sc->energy_up[ii][1];

              if (sc->f)
                en += sc->f(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
            }

            if (fij == en) {
              *i                            = ii + 1;
              *j                            = u;
              *k                            = u + 2;
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        } else {
          if (evaluate(ii, length, ii + 1, u, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            mm3   = (u < length) ? S1[u + 1] : -1;
            type  = get_pair_type_window(ii + 1, u, ptype);

            en = c[ii + 1][u - ii - 1] +
                 E_ExtLoop(type, mm5, mm3, P);

            if (sc) {
              if (sc->energy_up) {
                en += sc->energy_up[ii][1];
                if (u < length)
                  en += sc->energy_up[u + 1][1];
              }

              if (sc->f)
                en += sc->f(ii, length, ii + 1, u, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            if (fij == en) {
              *i                            = ii + 1;
              *j                            = u;
              *k                            = u + 2;
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
          type = get_pair_type_window(ii + 1, u, ptype);

          en = c[ii + 1][u - ii - 1] +
               E_ExtLoop(type, mm5, -1, P) +
               f3[u + 1];

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[ii][1];

            if (sc->f)
              en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
          }

          if (fij == en) {
            *i                            = ii + 1;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii + 1;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }

        if (u + 2 <= length) {
          if (evaluate(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3   = S1[u + 1];
            type  = get_pair_type_window(ii, u, ptype);

            en = c[ii][u - ii] +
                 E_ExtLoop(type, -1, mm3, P) +
                 f3[u + 2];

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[u + 1][1];

              if (sc->f)
                en += sc->f(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
            }

            if (fij == en) {
              *i                            = ii;
              *j                            = u;
              *k                            = u + 2;
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        } else {
          if (evaluate(ii, length, ii, u, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            mm3   = (u < length) ? S1[u + 1] : -1;
            type  = get_pair_type_window(ii, u, ptype);

            en = c[ii][u - ii] +
                 E_ExtLoop(type, -1, mm3, P);

            if (sc) {
              if ((sc->energy_up) && (u < length))
                en += sc->energy_up[u + 1][1];

              if (sc->f)
                en += sc->f(ii, length, ii, u, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            if (fij == en) {
              *i                            = ii;
              *j                            = u;
              *k                            = u + 2;
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          type = get_pair_type_window(ii, u, ptype);

          en = c[ii][u - ii] +
               E_ExtLoop(type, -1, -1, P) + f3[u + 1];

          if (sc)
            if (sc->f)
              en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }

      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f3_comparative(vrna_fold_compound_t *vc,
                           int                  *k,
                           int                  maxdist,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count)
{
  short                     **S, **S5, **S3;
  int                       length, fij, cc, fj, ii, u, *f3, **c, **ggg,
                            dangle_model, turn, with_gquad, ss, n_seq, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = vc->length;
  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3            = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  f3            = vc->matrices->f3_local;
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;

  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.cp         = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  ii = *k;

  /* nibble off unpaired 5' bases */
  do {
    fij = f3[ii];
    fj  = INF;

    if (evaluate(ii, maxdist, ii + 1, maxdist, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fj = f3[ii + 1];
      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fj += scs[ss]->energy_up[ii][1];

            if (scs[ss]->f)
              fj += scs[ss]->f(ii, maxdist, ii + 1, maxdist, VRNA_DECOMP_EXT_EXT, scs[ss]->data);
          }
      }
    }

    if (++ii > maxdist)
      break;
  } while (fij == fj);
  ii--;

  if (ii > maxdist - turn + 1) {
    /* no more pairs */
    *i  = *j = -1;
    *k  = 0;
    return 1;
  }

  /*
   *  must have found a decomposition
   *  i is paired. Find pairing partner
   */
  switch (dangle_model) {
    /* no dangles */
    case 0:
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            vrna_BT_gquad_mfe(vc, ii, u, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(ii, maxdist, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          cc = c[ii][u - ii];
          for (ss = 0; ss < n_seq; ss++) {
            type  = get_pair_type_md(S[ss][ii], S[ss][u], md);
            cc    += E_ExtLoop(type, -1, -1, P);
          }

          if (fij == cc + f3[u + 1]) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            vrna_BT_gquad_mfe(vc, ii, u, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(ii, maxdist, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          cc = c[ii][u - ii];
          for (ss = 0; ss < n_seq; ss++) {
            type  = get_pair_type_md(S[ss][ii], S[ss][u], md);
            cc    += E_ExtLoop(type, (ii > 1) ? S5[ss][ii] : -1, (u < length) ? S3[ss][u] : -1, P);
          }

          if (fij == cc + f3[u + 1]) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  int                   *i,
                  int                   maxj)
{
  int j, start;

  j     = -1;
  start = *i;

  if (fc) {
    char                      **ptype;
    short                     *S1;
    int                       traced2, length, turn, dangle_model, with_gquad, maxdist, type, cc,
                              **c, **ggg, *f3, fij, ii;
    vrna_param_t              *P;
    vrna_md_t                 *md;
    vrna_hc_t                 *hc;
    vrna_sc_t                 *sc;
    vrna_callback_hc_evaluate *evaluate;
    struct default_data       hc_dat_local;

    length        = fc->length;
    S1            = fc->sequence_encoding;
    ptype         = fc->ptype_local;
    f3            = fc->matrices->f3_local;
    c             = fc->matrices->c_local;
    ggg           = fc->matrices->ggg_local;
    hc            = fc->hc;
    sc            = fc->sc;
    P             = fc->params;
    md            = &(P->model_details);
    turn          = md->min_loop_size;
    dangle_model  = md->dangles;
    with_gquad    = md->gquad;
    maxdist       = MIN2(fc->window_size, maxj);
    traced2       = 0;
    ii            = start;

    hc_dat_local.mx_window  = hc->matrix_local;
    hc_dat_local.hc_up      = hc->up_ext;
    hc_dat_local.cp         = fc->cutpoint;

    if (hc->f) {
      evaluate            = &hc_default_user_window;
      hc_dat_local.hc_f   = hc->f;
      hc_dat_local.hc_dat = hc->data;
    } else {
      evaluate = &hc_default_window;
    }

    fij = f3[start];

    /* try to nibble off unpaired 5' bases */
    if ((sc) && (evaluate(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local))) {
      cc = f3[start + 1];

      if (sc->energy_up)
        cc += sc->energy_up[start][1];

      if (sc->f)
        cc += sc->f(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);

      if (fij == cc)
        /* simple 5' unpaired extensions, so we skip this hit */
        return 0;
    }

    /* get pairing partner j */
    switch (dangle_model) {
      case 0:
        for (j = start + turn; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = get_pair_type_window(start, j, ptype);

            cc = c[start][j - start] +
                 E_ExtLoop(type, -1, -1, P) +
                 f3[j + 1];

            if ((sc) && (sc->f))
              cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] + f3[j + 1];
            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }
        break;

      case 2:
        for (j = start + turn; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = get_pair_type_window(start, j, ptype);

            cc = c[start][j - start] +
                 E_ExtLoop(type, (start > 1) ? S1[start - 1] : -1,
                           (j < length) ? S1[j + 1] : -1, P) +
                 f3[j + 1];

            if ((sc) && (sc->f))
              cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] +
                 f3[j + 1];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }
        break;

      default:
        for (j = start + turn; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = get_pair_type_window(start, j, ptype);

            cc = c[start][j - start] +
                 E_ExtLoop(type, -1, -1, P) +
                 f3[j + 1];

            if ((sc) && (sc->f))
              cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (j < length) {
            if (evaluate(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
              cc = c[start][j - start] +
                   E_ExtLoop(type, -1, S1[j + 1], P) +
                   f3[j + 2];

              if (sc) {
                if (sc->energy_up)
                  cc += sc->energy_up[j + 1][1];

                if (sc->f)
                  cc += sc->f(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
              }

              if (fij == cc) {
                traced2 = 1;
                break;
              }
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] +
                 f3[j + 1];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            type = get_pair_type_window(start + 1, j, ptype);

            cc = c[start + 1][j - (start + 1)] +
                 E_ExtLoop(type, S1[start], -1, P) +
                 f3[j + 1];

            if (sc) {
              if (sc->energy_up)
                cc += sc->energy_up[start][1];

              if (sc->f)
                cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (j < length) {
            if (evaluate(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
              cc = c[start + 1][j - (start + 1)] +
                   E_ExtLoop(type, S1[start], S1[j + 1], P) +
                   f3[j + 2];

              if (sc) {
                if (sc->energy_up)
                  cc += sc->energy_up[start][1] +
                        sc->energy_up[j + 1][1];

                if (sc->f)
                  cc += sc->f(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
              }

              if (fij == cc) {
                traced2 = 1;
                break;
              }
            }
          }
        }
        break;
    }

    if (!traced2)
      j = -1;
  }

  *i = start;
  return j;
}


PRIVATE int
BT_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              int                   *i,
                              int                   maxj)
{
  int j, start;

  j = -1;

  if (fc) {
    short                     **S, **S5, **S3;
    int                       traced2, length, turn, dangle_model, with_gquad, maxdist, cc, **c,
                              **ggg, *f3, s, tt, n_seq, fij;
    vrna_param_t              *P;
    vrna_md_t                 *md;
    vrna_hc_t                 *hc;
    vrna_sc_t                 **scs;
    vrna_callback_hc_evaluate *evaluate;
    struct default_data       hc_dat_local;

    length        = fc->length;
    n_seq         = fc->n_seq;
    S             = fc->S;
    S5            = fc->S5;   /* S5[s][start] holds next base 5' of start in sequence s */
    S3            = fc->S3;   /* Sl[s][start] holds next base 3' of start in sequence s */
    f3            = fc->matrices->f3_local;
    c             = fc->matrices->c_local;
    ggg           = fc->matrices->ggg_local;
    hc            = fc->hc;
    scs           = fc->scs;
    P             = fc->params;
    md            = &(P->model_details);
    turn          = md->min_loop_size;
    dangle_model  = md->dangles;
    with_gquad    = md->gquad;
    maxdist       = MIN2(fc->window_size, maxj);
    traced2       = 0;
    start         = *i;

    hc_dat_local.mx_window  = hc->matrix_local;
    hc_dat_local.hc_up      = hc->up_ext;
    hc_dat_local.cp         = fc->cutpoint;

    if (hc->f) {
      evaluate            = &hc_default_user_window;
      hc_dat_local.hc_f   = hc->f;
      hc_dat_local.hc_dat = hc->data;
    } else {
      evaluate = &hc_default_window;
    }

    fij = f3[start];

    /* try to nibble off unpaired 5' bases */
    if ((scs) && (evaluate(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local))) {
      cc = f3[start + 1];

      for (s = 0; s < n_seq; s++)
        if (scs[s]) {
          if (scs[s]->energy_up)
            cc += scs[s]->energy_up[start][1];

          if (scs[s]->f)
            cc +=
              scs[s]->f(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, scs[s]->data);
        }

      if (fij == cc)
        /* simple 5' unpaired extensions, so we skip this hit */
        return 0;
    }

    /* get pairing partner j */
    switch (dangle_model) {
      case 0:
        for (j = start + turn; j <= MIN2(start + maxdist, length - 1); j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            cc = c[start][j - start];

            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][start], S[s][j], md);
              cc  += E_ExtLoop(tt, -1, -1, P);
            }

            if (fij == cc + f3[j + 1]) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc + f3[j + 1]) {
              traced2 = 1;
              break;
            }
          }
        }

        if ((!traced2) && (length <= start + maxdist)) {
          j = length;
          if (evaluate(start, length, start, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            cc = c[start][j - start];

            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][start], S[s][j], md);
              cc  += E_ExtLoop(tt, -1, -1, P);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }

        break;

      case 2:
        for (j = start + turn; j <= MIN2(start + maxdist, length - 1); j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            cc = c[start][j - start];

            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][start], S[s][j], md);
              cc  += E_ExtLoop(tt, (start > 1) ? S5[s][start] : -1, S3[s][j], P);
            }

            if (fij == cc + f3[j + 1]) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc + f3[j + 1]) {
              traced2 = 1;
              break;
            }
          }
        }

        if ((!traced2) && (length <= start + maxdist)) {
          j = length;
          if (evaluate(start, length, start, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            cc = c[start][j - start];
            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][start], S[s][j], md);
              cc  += E_ExtLoop(tt, (start > 1) ? S5[s][start] :  -1, -1, P);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }

        break;
    }

    if (!traced2)
      j = -1;
  }

  return j;
}


PUBLIC vrna_mx_pf_aux_el_t *
vrna_exp_E_ext_fast_init(vrna_fold_compound_t *vc)
{
  vrna_mx_pf_aux_el_t *aux_mx = NULL;

  if (vc) {
    unsigned int              u, s;
    int                       i, j, max_j, d, n, turn, ij, *idx, *iidx, *hc_up;
    FLT_OR_DBL                *q, **q_local, *scale;
    vrna_callback_hc_evaluate *evaluate;
    struct default_data       hc_dat_local;

    n     = (int)vc->length;
    idx   = vc->jindx;
    iidx  = vc->iindx;
    turn  = vc->exp_params->model_details.min_loop_size;
    scale = vc->exp_matrices->scale;
    hc_up = vc->hc->up_ext;

    if (vc->hc->type == VRNA_HC_WINDOW) {
      hc_dat_local.mx_window  = vc->hc->matrix_local;
      hc_dat_local.hc_up      = hc_up;
      hc_dat_local.cp         = vc->cutpoint;

      if (vc->hc->f) {
        evaluate            = &hc_default_user_window;
        hc_dat_local.hc_f   = vc->hc->f;
        hc_dat_local.hc_dat = vc->hc->data;
      } else {
        evaluate = &hc_default_window;
      }
    } else {
      hc_dat_local.idx    = idx;
      hc_dat_local.mx     = vc->hc->matrix;
      hc_dat_local.hc_up  = hc_up;
      hc_dat_local.cp     = vc->cutpoint;

      if (vc->hc->f) {
        evaluate            = &hc_default_user;
        hc_dat_local.hc_f   = vc->hc->f;
        hc_dat_local.hc_dat = vc->hc->data;
      } else {
        evaluate = &hc_default;
      }
    }

    /* allocate memory for helper arrays */
    aux_mx            = (vrna_mx_pf_aux_el_t *)vrna_alloc(sizeof(vrna_mx_pf_aux_el_t));
    aux_mx->qq        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qq1       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqu_size  = 0;
    aux_mx->qqu       = NULL;

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      vrna_sc_t *sc         = vc->sc;
      vrna_ud_t *domains_up = vc->domains_up;
      int       with_ud     = (domains_up && domains_up->exp_energy_cb);

      /* pre-processing ligand binding production rule(s) and auxiliary memory */
      if (with_ud) {
        int ud_max_size = 0;
        for (u = 0; u < domains_up->uniq_motif_count; u++)
          if (ud_max_size < domains_up->uniq_motif_size[u])
            ud_max_size = domains_up->uniq_motif_size[u];

        aux_mx->qqu_size  = ud_max_size;
        aux_mx->qqu       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));

        for (u = 0; u <= ud_max_size; u++)
          aux_mx->qqu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
      }

      if (vc->hc->type == VRNA_HC_WINDOW) {
        q_local = vc->exp_matrices->q_local;
        max_j   = MIN2(turn + 1, vc->window_size);
        max_j   = MIN2(max_j, n);
        for (j = 1; j <= max_j; j++)
          for (i = 1; i <= j; i++) {
            if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
              q_local[i][j] = scale[(j - i + 1)];
              if (sc) {
                if (sc->exp_energy_up)
                  q_local[i][j] *= sc->exp_energy_up[i][j - i + 1];

                if (sc->exp_f)
                  q_local[i][j] *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
              }

              if (with_ud) {
                q_local[i][j] += q_local[i][j] *
                                 domains_up->exp_energy_cb(vc,
                                                           i, j,
                                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                                           domains_up->data);
              }
            } else {
              q_local[i][j] = 0;
            }
          }
      } else {
        q = vc->exp_matrices->q;
        for (d = 0; d <= turn; d++)
          for (i = 1; i <= n - d; i++) {
            j   = i + d;
            ij  = iidx[i] - j;

            if (j > n)
              continue;

            if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
              q[ij] = scale[d + 1];

              if (sc) {
                if (sc->exp_energy_up)
                  q[ij] *= sc->exp_energy_up[i][d + 1];

                if (sc->exp_f)
                  q[ij] *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
              }

              if (with_ud) {
                q[ij] += q[ij] *
                         domains_up->exp_energy_cb(vc,
                                                   i, j,
                                                   VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                                   domains_up->data);
              }
            } else {
              q[ij] = 0.;
            }
          }
      }
    } else if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
      vrna_sc_t       **scs = vc->scs;
      unsigned int    **a2s = vc->a2s;
      q = vc->exp_matrices->q;
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;
          if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
            q[ij] = scale[d + 1];

            if (scs) {
              for (s = 0; s < vc->n_seq; s++)
                if (scs[s]) {
                  u = d + 1 /* a2s[s][j] - a2s[s][i] + 1 */;
                  if (scs[s]->exp_energy_up)
                    q[ij] *= scs[s]->exp_energy_up[a2s[s][i]][u];
                }
            }
          } else {
            q[ij] = 0.;
          }
        }
    }
  }

  return aux_mx;
}


PUBLIC void
vrna_exp_E_ext_fast_rotate(vrna_fold_compound_t *vc,
                           vrna_mx_pf_aux_el_t  *aux_mx)
{
  if (vc && aux_mx) {
    int         u;
    FLT_OR_DBL  *tmp;

    tmp         = aux_mx->qq1;
    aux_mx->qq1 = aux_mx->qq;
    aux_mx->qq  = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqu) {
      tmp = aux_mx->qqu[aux_mx->qqu_size];
      for (u = aux_mx->qqu_size; u > 0; u--)
        aux_mx->qqu[u] = aux_mx->qqu[u - 1];
      aux_mx->qqu[0] = tmp;
    }
  }
}


PUBLIC void
vrna_exp_E_ext_fast_free(vrna_fold_compound_t *vc,
                         vrna_mx_pf_aux_el_t  *aux_mx)
{
  if (vc && aux_mx) {
    int u;

    free(aux_mx->qq);
    free(aux_mx->qq1);

    if (aux_mx->qqu) {
      for (u = 0; u <= aux_mx->qqu_size; u++)
        free(aux_mx->qqu[u]);

      free(aux_mx->qqu);
    }

    free(aux_mx);
  }
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ext_fast(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    vrna_mx_pf_aux_el_t   *aux_mx)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return exp_E_ext_fast_window(vc, i, j, aux_mx);
        else
          return exp_E_ext_fast(vc, i, j, aux_mx);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return exp_E_ext_fast_comparative(vc, i, j, aux_mx);
        break;

      default:
        vrna_message_warning("vrna_exp_E_ext_fast@exterior_loops.c: Unknown fold_compound type");
        return 0.;
        break;
    }
  } else {
    return 0.;
  }
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_el_t  *aux_mx)
{
  short                     *S1, *S2;
  int                       n, *iidx, k, ij, kl, with_ud, u, circular, with_gquad, type;
  FLT_OR_DBL                qbt1, *q, *qb, *qq, *qq1, **qqu, q_temp, *scale, q_temp2, *G;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n                   = (int)vc->length;
  iidx                = vc->iindx;
  ij                  = iidx[i] - j;
  qq                  = aux_mx->qq;
  qq1                 = aux_mx->qq1;
  qqu                 = aux_mx->qqu;
  q                   = vc->exp_matrices->q;
  qb                  = vc->exp_matrices->qb;
  G                   = vc->exp_matrices->G;
  scale               = vc->exp_matrices->scale;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  sc                  = vc->sc;
  domains_up          = vc->domains_up;
  circular            = md->circ;
  with_gquad          = md->gquad;
  with_ud             = (domains_up && domains_up->exp_energy_cb);
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            q_temp2 = qqu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      scale[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    S1      = vc->sequence_encoding;
    S2      = vc->sequence_encoding2;
    type    = get_pair_type_md(S2[i], S2[j], md);
    q_temp  = qb[ij] *
              exp_E_ExtLoop(type,
                            ((i > 1) || circular) ? S1[i - 1] : -1,
                            ((j < n) || circular) ? S1[j + 1] : -1,
                            pf_params);

    if (sc)
      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

    qbt1 += q_temp;
  }

  if (with_gquad)
    qbt1 += G[ij];

  qq[i] = qbt1;

  if (with_ud)
    qqu[0][i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[i][u];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
    }

    qbt1 += q_temp;

    if (with_ud) {
      qbt1 += q_temp *
              domains_up->exp_energy_cb(vc,
                                        i, j,
                                        VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                        domains_up->data);
    }
  }

  kl = iidx[i] - j + 1;
  if (sc && sc->exp_f) {
    for (k = j; k > i; k--, kl++) {
      q_temp = q[kl] *
               qq[k] *
               sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, sc->data);

      qbt1 += q_temp;
    }
  } else {
    for (k = j; k > i; k--, kl++)
      qbt1 += q[kl] *
              qq[k];
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      vrna_mx_pf_aux_el_t   *aux_mx)
{
  short                     *S1;
  int                       n, k, with_ud, u, circular, with_gquad, type, winSize, turn;
  FLT_OR_DBL                qbt1, **q, **qb, *qq, *qq1, **qqu, q_temp, *scale, q_temp2, **G;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n           = (int)vc->length;
  winSize     = vc->window_size;
  qq          = aux_mx->qq;
  qq1         = aux_mx->qq1;
  qqu         = aux_mx->qqu;
  q           = vc->exp_matrices->q_local;
  qb          = vc->exp_matrices->qb_local;
  G           = vc->exp_matrices->G_local;
  scale       = vc->exp_matrices->scale;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  domains_up  = vc->domains_up;
  circular    = md->circ;
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  with_ud     = (domains_up && domains_up->exp_energy_cb);

  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.cp         = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /*
   *  init exterior loop contributions for small segments [i, j]
   *  that can only be unpaired.
   *  We do this only once for the very first segment ending at j
   */
  if (i == j - turn - 1) {
    for (k = j; k >= MAX2(1, j - turn); k--) {
      if (evaluate(k, j, k, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
        q[k][j] = scale[(j - k + 1)];
        if (sc) {
          if (sc->exp_energy_up)
            q[k][j] *= sc->exp_energy_up[k][j - k + 1];

          if (sc->exp_f)
            q[k][j] *= sc->exp_f(k, j, k, j, VRNA_DECOMP_EXT_UP, sc->data);
        }

        if (with_ud) {
          q[k][j] += q[k][j] *
                     domains_up->exp_energy_cb(vc,
                                               k, j,
                                               VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                               domains_up->data);
        }
      } else {
        q[k][j] = 0;
      }
    }
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            q_temp2 = qqu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      scale[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    S1      = vc->sequence_encoding;
    type    = get_pair_type_md(S1[i], S1[j], md);
    q_temp  = qb[i][j] *
              exp_E_ExtLoop(type,
                            ((i > 1) || circular) ? S1[i - 1] : -1,
                            ((j < n) || circular) ? S1[j + 1] : -1,
                            pf_params);

    if (sc)
      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

    qbt1 += q_temp;
  }

  if (with_gquad)
    qbt1 += G[i][j];

  qq[i] = qbt1;

  if (with_ud)
    qqu[0][i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[i][u];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
    }

    qbt1 += q_temp;

    if (with_ud) {
      qbt1 += q_temp *
              domains_up->exp_energy_cb(vc,
                                        i, j,
                                        VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                        domains_up->data);
    }
  }

  if (sc && sc->exp_f) {
    for (k = j; k > i; k--) {
      q_temp = q[i][k - 1] *
               qq[k] *
               sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, sc->data);

      qbt1 += q_temp;
    }
  } else {
    for (k = j; k > i; k--)
      qbt1 += q[i][k - 1] *
              qq[k];
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           vrna_mx_pf_aux_el_t  *aux_mx)
{
  int                       n, s, n_seq, *iidx, k, ij, kl, u, circular, type;
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  FLT_OR_DBL                qbt1, *q, *qb, *qq, *qq1, q_temp, *scale;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n                   = (int)vc->length;
  n_seq               = vc->n_seq;
  iidx                = vc->iindx;
  ij                  = iidx[i] - j;
  S                   = vc->S;
  S5                  = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
  S3                  = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s                 = vc->a2s;
  qq                  = aux_mx->qq;
  qq1                 = aux_mx->qq1;
  q                   = vc->exp_matrices->q;
  qb                  = vc->exp_matrices->qb;
  scale               = vc->exp_matrices->scale;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  scs                 = vc->scs;
  circular            = md->circ;
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] *
             scale[1];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][j]][1];
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    q_temp = qb[ij];

    for (s = 0; s < n_seq; s++) {
      type    = get_pair_type_md(S[s][i], S[s][j], md);
      q_temp  *= exp_E_ExtLoop(type,
                               ((i > 1) || circular) ? S5[s][i] : -1,
                               ((j < n) || circular) ? S3[s][j] : -1,
                               pf_params);
    }

    qbt1 += q_temp;
  }

  qq[i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][i]][a2s[s][j] - a2s[s][i] + 1];
      }
    }

    qbt1 += q_temp;
  }

  kl = iidx[i] - j + 1;
  for (k = j; k > i; k--, kl++)
    qbt1 += q[kl] *
            qq[k];

  return qbt1;
}


PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data)
{
  int                 kl, di, dj;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      kl = dat->idx[j] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      kl = dat->idx[j - 1] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;
          if (dat->hc_up[j] == 0)
            eval = (unsigned char)0;

          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    default:
      vrna_message_error("hc_cb@exterior_loops.c: Unrecognized decomposition %d", d);
  }
  return eval;
}


PRIVATE unsigned char
hc_default_window(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data)
{
  int                 di, dj;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      if (dat->mx_window[l][j - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      if (dat->mx_window[i][k - i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (j != k) {
          /* otherwise, stem spans from i to j */
          dj = l - k - 1;
          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (dat->mx_window[l][j - 1 - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (dat->mx_window[i + 1][k - (i + 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;

          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    default:
      vrna_message_error("hc_cb@exterior_loops.c: Unrecognized decomposition %d", d);
  }
  return eval;
}


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_default_user_window(int            i,
                       int            j,
                       int            k,
                       int            l,
                       unsigned char  d,
                       void           *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default_window(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}
