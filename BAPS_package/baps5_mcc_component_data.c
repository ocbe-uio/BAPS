/*
 * MATLAB Compiler: 4.13 (R2010a)
 * Date: Fri Jun  8 14:40:22 2012
 * Arguments: "-B" "macro_default" "-m" "-W" "main" "-T" "link:exe"
 * "./general/baps5.m" "-a" "./admixture" "-a" "./general" "-a" "./graph" "-a"
 * "./independent" "-a" "./linkage" "-a" "./parallel" "-a" "./spatial" "-d"
 * "./BAPS_package" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_baps5_session_key[] = {
    '9', '1', '9', '3', '4', '2', '1', 'F', 'D', 'E', '6', '3', '3', '4', 'C',
    '6', 'B', '4', '9', '1', 'B', 'D', '5', 'A', '9', '5', 'A', '1', '1', '8',
    '1', '5', '3', 'F', '9', 'A', '9', 'D', '4', 'B', '7', 'F', 'E', '6', '4',
    'D', 'D', '4', 'E', '4', 'D', '7', 'B', 'F', '4', '3', 'F', 'C', 'D', '1',
    '7', '0', '5', '2', '1', '4', 'F', '9', 'C', '3', '4', '8', '6', 'F', '5',
    'D', 'A', 'A', 'B', 'D', '4', '3', '1', '3', '9', '2', '1', 'F', 'C', '1',
    'F', 'B', 'F', 'B', 'F', '6', 'C', '2', '7', '9', 'E', 'A', '7', '0', 'E',
    '2', '0', 'D', '7', 'E', 'E', 'A', '6', 'F', '1', '1', '7', '0', 'F', '0',
    '7', '9', '7', '4', '4', 'F', 'B', 'B', 'F', '5', 'D', '7', 'D', 'D', '0',
    '8', 'B', '9', '3', '4', '6', '5', 'F', '1', '9', '1', 'B', '9', '1', '8',
    '3', 'C', 'B', '6', 'C', '9', '7', '5', '4', 'D', '5', '4', 'C', 'F', '6',
    'C', 'D', '3', 'A', 'C', '0', '6', '3', '6', '3', 'C', '5', '6', 'A', '9',
    'C', 'D', 'D', '7', 'A', 'E', 'F', '0', '7', 'E', 'C', 'C', '3', '1', '0',
    '1', 'C', 'C', '1', '0', '9', '6', 'C', '0', '2', '3', '0', '2', '8', 'A',
    '7', 'D', 'A', '3', '0', '9', '0', '0', '4', '5', '9', '2', '8', '5', 'D',
    'B', '5', '4', '2', '2', '0', '4', 'B', '7', '0', '0', 'A', 'F', '4', '8',
    '2', '7', 'D', '5', '5', '0', '9', '3', '8', 'A', '5', 'D', 'F', '7', '1',
    'E', '\0'};

const unsigned char __MCC_baps5_public_key[] = {
    '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9', '2',
    'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1', '0', '1',
    '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B', '0', '0', '3',
    '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1', '0', '0', 'C', '4',
    '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3', 'A', '5', '2', '0', '6',
    '5', '8', 'F', '6', 'F', '8', 'E', '0', '1', '3', '8', 'C', '4', '3', '1',
    '5', 'B', '4', '3', '1', '5', '2', '7', '7', 'E', 'D', '3', 'F', '7', 'D',
    'A', 'E', '5', '3', '0', '9', '9', 'D', 'B', '0', '8', 'E', 'E', '5', '8',
    '9', 'F', '8', '0', '4', 'D', '4', 'B', '9', '8', '1', '3', '2', '6', 'A',
    '5', '2', 'C', 'C', 'E', '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4',
    'D', '0', '8', '5', 'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2',
    'E', 'D', 'E', '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6',
    '3', '7', '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E',
    '6', '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
    '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1', 'B',
    'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9', '9', '0',
    '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0', 'B', '6', '1',
    'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B', '5', '8', 'F', 'C',
    '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6', 'E', 'B', '7', 'E', 'C',
    'D', '3', '1', '7', '8', 'B', '5', '6', 'A', 'B', '0', 'F', 'A', '0', '6',
    'D', 'D', '6', '4', '9', '6', '7', 'C', 'B', '1', '4', '9', 'E', '5', '0',
    '2', '0', '1', '1', '1', '\0'};

static const char * MCC_baps5_matlabpath_data[] = 
  { "baps5/", "$TOOLBOXDEPLOYDIR/", "general/", "admixture/",
    "graph/", "independent/", "linkage/", "parallel/", "spatial/",
    "$TOOLBOXMATLABDIR/general/", "$TOOLBOXMATLABDIR/ops/",
    "$TOOLBOXMATLABDIR/lang/", "$TOOLBOXMATLABDIR/elmat/",
    "$TOOLBOXMATLABDIR/randfun/", "$TOOLBOXMATLABDIR/elfun/",
    "$TOOLBOXMATLABDIR/specfun/", "$TOOLBOXMATLABDIR/matfun/",
    "$TOOLBOXMATLABDIR/datafun/", "$TOOLBOXMATLABDIR/polyfun/",
    "$TOOLBOXMATLABDIR/funfun/", "$TOOLBOXMATLABDIR/sparfun/",
    "$TOOLBOXMATLABDIR/scribe/", "$TOOLBOXMATLABDIR/graph2d/",
    "$TOOLBOXMATLABDIR/graph3d/", "$TOOLBOXMATLABDIR/specgraph/",
    "$TOOLBOXMATLABDIR/graphics/", "$TOOLBOXMATLABDIR/uitools/",
    "$TOOLBOXMATLABDIR/strfun/", "$TOOLBOXMATLABDIR/imagesci/",
    "$TOOLBOXMATLABDIR/iofun/", "$TOOLBOXMATLABDIR/audiovideo/",
    "$TOOLBOXMATLABDIR/timefun/", "$TOOLBOXMATLABDIR/datatypes/",
    "$TOOLBOXMATLABDIR/verctrl/", "$TOOLBOXMATLABDIR/codetools/",
    "$TOOLBOXMATLABDIR/helptools/", "$TOOLBOXMATLABDIR/demos/",
    "$TOOLBOXMATLABDIR/timeseries/", "$TOOLBOXMATLABDIR/hds/",
    "$TOOLBOXMATLABDIR/guide/", "$TOOLBOXMATLABDIR/plottools/",
    "toolbox/local/", "$TOOLBOXMATLABDIR/datamanager/",
    "toolbox/bioinfo/bioinfo/", "toolbox/bioinfo/biolearning/",
    "toolbox/compiler/", "toolbox/stats/" };

static const char * MCC_baps5_classpath_data[] = 
  { "" };

static const char * MCC_baps5_libpath_data[] = 
  { "" };

static const char * MCC_baps5_app_opts_data[] = 
  { "" };

static const char * MCC_baps5_run_opts_data[] = 
  { "" };

static const char * MCC_baps5_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_baps5_component_data = { 

  /* Public key data */
  __MCC_baps5_public_key,

  /* Component name */
  "baps5",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_baps5_session_key,

  /* Component's MATLAB Path */
  MCC_baps5_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  47,

  /* Component's Java class path */
  MCC_baps5_classpath_data,
  /* Number of directories in the Java class path */
  0,

  /* Component's load library path (for extra shared libraries) */
  MCC_baps5_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_baps5_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_baps5_run_opts_data,
  /* Number of MCR global runtime options */
  0,
  
  /* Component preferences directory */
  "baps5_966BADC4A6C056D5FFF68D4F39D76F83",

  /* MCR warning status data */
  MCC_baps5_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


