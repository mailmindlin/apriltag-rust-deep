#include "apriltag-src/apriltag.h"
#include "apriltag-src/apriltag_math.h"
#include "apriltag-src/common/unionfind.h"
#include "apriltag-src/common/zarray.h"
#include "apriltag-src/common/zhash.h"
#include "apriltag-src/common/zmaxheap.h"
#include "apriltag-src/common/homography.h"

extern zarray_t *apriltag_quad_thresh(apriltag_detector_t *td, image_u8_t *im);
unionfind_t* connected_components(apriltag_detector_t *td, image_u8_t* threshim, int w, int h, int ts);
zarray_t* do_gradient_clusters(image_u8_t* threshim, int ts, int y0, int y1, int w, int nclustermap, unionfind_t* uf, zarray_t* clusters);
zarray_t* merge_clusters(zarray_t* c1, zarray_t* c2);
zarray_t* gradient_clusters(apriltag_detector_t *td, image_u8_t* threshim, int w, int h, int ts, unionfind_t* uf);
zarray_t* fit_quads(apriltag_detector_t *td, int w, int h, zarray_t* clusters, image_u8_t* im);
struct pt
{
    // Note: these represent 2*actual value.
    uint16_t x, y;
    int16_t gx, gy;

    float slope;
};

// return 1 if the quad looks okay, 0 if it should be discarded
int fit_quad(
        apriltag_detector_t *td,
        image_u8_t *im,
        zarray_t *cluster,
        struct quad *quad,
        int tag_width,
        bool normal_border,
        bool reversed_border);

struct line_fit_pt
{
    double Mx, My;
    double Mxx, Myy, Mxy;
    double W; // total weight
};
void fit_line(struct line_fit_pt *lfps, int sz, int i0, int i1, double *lineparm, double *err, double *mse);
struct line_fit_pt* compute_lfps(int sz, zarray_t* cluster, image_u8_t* im);

void refine_edges(apriltag_detector_t *td, image_u8_t *im_orig, struct quad *quad);
int quad_update_homographies(struct quad *quad);

struct quick_decode_entry
{
    uint64_t rcode;   // the queried code
    uint16_t id;      // the tag ID (a small integer)
    uint8_t hamming;  // how many errors corrected?
    uint8_t rotation; // number of rotations [0, 3]
};
float quad_decode(apriltag_detector_t* td, apriltag_family_t *family, image_u8_t *im, struct quad *quad, struct quick_decode_entry *entry, image_u8_t *im_samples);

void quick_decode_init(apriltag_family_t *family, int maxhamming);
void quick_decode_uninit(apriltag_family_t *fam);

int quad_segment_maxima(apriltag_detector_t *td, zarray_t *cluster, struct line_fit_pt *lfps, int indices[4]);
void ptsort(struct pt *pts, int sz);