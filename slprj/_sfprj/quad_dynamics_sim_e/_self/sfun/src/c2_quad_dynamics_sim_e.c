/* Include files */

#include <stddef.h>
#include "blas.h"
#include "quad_dynamics_sim_e_sfun.h"
#include "c2_quad_dynamics_sim_e.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "quad_dynamics_sim_e_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c2_debug_family_names[61] = { "m", "g", "l", "Jxx", "Jyy",
  "Jzz", "k1", "k2", "k3", "k4", "b1", "b2", "b3", "b4", "w1", "w2", "w3", "w4",
  "x", "x_dot", "y", "y_dot", "z", "z_dot", "phi", "phi_dot", "theta",
  "theta_dot", "psi", "psi_dot", "multipliers", "upphi_x", "upphi_y", "upphi_z",
  "T_e2E", "force_mult_1", "force_mult_2", "dx1_dphi", "dx2_dphi", "upphi_phi",
  "dx1_dtheta", "dx2_dtheta", "upphi_theta", "dx1_dpsi", "dx2_dpsi", "upphi_psi",
  "delta1", "delta2", "delta3", "ddot", "nargin", "nargout", "states", "params",
  "w", "x_ddot", "y_ddot", "z_ddot", "phi_ddot", "theta_ddot", "psi_ddot" };

/* Function Declarations */
static void initialize_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void initialize_params_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void enable_c2_quad_dynamics_sim_e(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance);
static void disable_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void c2_update_debugger_state_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void set_sim_state_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance, const mxArray *c2_st);
static void finalize_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void sf_c2_quad_dynamics_sim_e(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance);
static void c2_chartstep_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void initSimStructsc2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void registerMessagesc2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static real_T c2_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_psi_ddot, const char_T *c2_identifier);
static real_T c2_b_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_c_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[3]);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_d_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[9]);
static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_info_helper(c2_ResolvedFunctionInfo c2_info[63]);
static real_T c2_mpower(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance,
  real_T c2_a);
static void c2_eml_scalar_eg(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance);
static void c2_b_eml_scalar_eg(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance);
static real_T c2_dot(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance,
                     real_T c2_a[3], real_T c2_b[3]);
static void c2_c_eml_scalar_eg(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance);
static void c2_mldivide(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance,
  real_T c2_A[9], real_T c2_B[3], real_T c2_Y[3]);
static void c2_eml_warning(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance);
static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_e_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_f_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_quad_dynamics_sim_e, const
  char_T *c2_identifier);
static uint8_T c2_g_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void init_dsm_address_info(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
  chartInstance->c2_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c2_is_active_c2_quad_dynamics_sim_e = 0U;
}

static void initialize_params_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
}

static void enable_c2_quad_dynamics_sim_e(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c2_update_debugger_state_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  real_T c2_hoistedGlobal;
  real_T c2_u;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_hoistedGlobal;
  real_T c2_b_u;
  const mxArray *c2_c_y = NULL;
  real_T c2_c_hoistedGlobal;
  real_T c2_c_u;
  const mxArray *c2_d_y = NULL;
  real_T c2_d_hoistedGlobal;
  real_T c2_d_u;
  const mxArray *c2_e_y = NULL;
  real_T c2_e_hoistedGlobal;
  real_T c2_e_u;
  const mxArray *c2_f_y = NULL;
  real_T c2_f_hoistedGlobal;
  real_T c2_f_u;
  const mxArray *c2_g_y = NULL;
  uint8_T c2_g_hoistedGlobal;
  uint8_T c2_g_u;
  const mxArray *c2_h_y = NULL;
  real_T *c2_phi_ddot;
  real_T *c2_psi_ddot;
  real_T *c2_theta_ddot;
  real_T *c2_x_ddot;
  real_T *c2_y_ddot;
  real_T *c2_z_ddot;
  c2_psi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c2_theta_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c2_phi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_z_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_y_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_x_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_st = NULL;
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellarray(7), FALSE);
  c2_hoistedGlobal = *c2_phi_ddot;
  c2_u = c2_hoistedGlobal;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_b_hoistedGlobal = *c2_psi_ddot;
  c2_b_u = c2_b_hoistedGlobal;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  c2_c_hoistedGlobal = *c2_theta_ddot;
  c2_c_u = c2_c_hoistedGlobal;
  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", &c2_c_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 2, c2_d_y);
  c2_d_hoistedGlobal = *c2_x_ddot;
  c2_d_u = c2_d_hoistedGlobal;
  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", &c2_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 3, c2_e_y);
  c2_e_hoistedGlobal = *c2_y_ddot;
  c2_e_u = c2_e_hoistedGlobal;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_e_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 4, c2_f_y);
  c2_f_hoistedGlobal = *c2_z_ddot;
  c2_f_u = c2_f_hoistedGlobal;
  c2_g_y = NULL;
  sf_mex_assign(&c2_g_y, sf_mex_create("y", &c2_f_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 5, c2_g_y);
  c2_g_hoistedGlobal = chartInstance->c2_is_active_c2_quad_dynamics_sim_e;
  c2_g_u = c2_g_hoistedGlobal;
  c2_h_y = NULL;
  sf_mex_assign(&c2_h_y, sf_mex_create("y", &c2_g_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 6, c2_h_y);
  sf_mex_assign(&c2_st, c2_y, FALSE);
  return c2_st;
}

static void set_sim_state_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance, const mxArray *c2_st)
{
  const mxArray *c2_u;
  real_T *c2_phi_ddot;
  real_T *c2_psi_ddot;
  real_T *c2_theta_ddot;
  real_T *c2_x_ddot;
  real_T *c2_y_ddot;
  real_T *c2_z_ddot;
  c2_psi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c2_theta_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c2_phi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_z_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_y_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_x_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c2_doneDoubleBufferReInit = TRUE;
  c2_u = sf_mex_dup(c2_st);
  *c2_phi_ddot = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c2_u, 0)), "phi_ddot");
  *c2_psi_ddot = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c2_u, 1)), "psi_ddot");
  *c2_theta_ddot = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c2_u, 2)), "theta_ddot");
  *c2_x_ddot = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u,
    3)), "x_ddot");
  *c2_y_ddot = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u,
    4)), "y_ddot");
  *c2_z_ddot = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u,
    5)), "z_ddot");
  chartInstance->c2_is_active_c2_quad_dynamics_sim_e = c2_f_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 6)),
     "is_active_c2_quad_dynamics_sim_e");
  sf_mex_destroy(&c2_u);
  c2_update_debugger_state_c2_quad_dynamics_sim_e(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
}

static void sf_c2_quad_dynamics_sim_e(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance)
{
  int32_T c2_i0;
  int32_T c2_i1;
  int32_T c2_i2;
  real_T *c2_x_ddot;
  real_T *c2_y_ddot;
  real_T *c2_z_ddot;
  real_T *c2_phi_ddot;
  real_T *c2_theta_ddot;
  real_T *c2_psi_ddot;
  real_T (*c2_w)[4];
  real_T (*c2_params)[14];
  real_T (*c2_states)[12];
  c2_w = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
  c2_params = (real_T (*)[14])ssGetInputPortSignal(chartInstance->S, 1);
  c2_psi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c2_theta_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c2_phi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_z_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_y_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_x_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_states = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c2_sfEvent);
  for (c2_i0 = 0; c2_i0 < 12; c2_i0++) {
    _SFD_DATA_RANGE_CHECK((*c2_states)[c2_i0], 0U);
  }

  _SFD_DATA_RANGE_CHECK(*c2_x_ddot, 1U);
  _SFD_DATA_RANGE_CHECK(*c2_y_ddot, 2U);
  _SFD_DATA_RANGE_CHECK(*c2_z_ddot, 3U);
  _SFD_DATA_RANGE_CHECK(*c2_phi_ddot, 4U);
  _SFD_DATA_RANGE_CHECK(*c2_theta_ddot, 5U);
  _SFD_DATA_RANGE_CHECK(*c2_psi_ddot, 6U);
  for (c2_i1 = 0; c2_i1 < 14; c2_i1++) {
    _SFD_DATA_RANGE_CHECK((*c2_params)[c2_i1], 7U);
  }

  for (c2_i2 = 0; c2_i2 < 4; c2_i2++) {
    _SFD_DATA_RANGE_CHECK((*c2_w)[c2_i2], 8U);
  }

  chartInstance->c2_sfEvent = CALL_EVENT;
  c2_chartstep_c2_quad_dynamics_sim_e(chartInstance);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_quad_dynamics_sim_eMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c2_chartstep_c2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
  int32_T c2_i3;
  real_T c2_states[12];
  int32_T c2_i4;
  real_T c2_params[14];
  int32_T c2_i5;
  real_T c2_w[4];
  uint32_T c2_debug_family_var_map[61];
  real_T c2_m;
  real_T c2_g;
  real_T c2_l;
  real_T c2_Jxx;
  real_T c2_Jyy;
  real_T c2_Jzz;
  real_T c2_k1;
  real_T c2_k2;
  real_T c2_k3;
  real_T c2_k4;
  real_T c2_b1;
  real_T c2_b2;
  real_T c2_b3;
  real_T c2_b4;
  real_T c2_w1;
  real_T c2_w2;
  real_T c2_w3;
  real_T c2_w4;
  real_T c2_x;
  real_T c2_x_dot;
  real_T c2_y;
  real_T c2_y_dot;
  real_T c2_z;
  real_T c2_z_dot;
  real_T c2_phi;
  real_T c2_phi_dot;
  real_T c2_theta;
  real_T c2_theta_dot;
  real_T c2_psi;
  real_T c2_psi_dot;
  real_T c2_multipliers[9];
  real_T c2_upphi_x;
  real_T c2_upphi_y;
  real_T c2_upphi_z;
  real_T c2_T_e2E[9];
  real_T c2_force_mult_1[3];
  real_T c2_force_mult_2[3];
  real_T c2_dx1_dphi[3];
  real_T c2_dx2_dphi[3];
  real_T c2_upphi_phi;
  real_T c2_dx1_dtheta[3];
  real_T c2_dx2_dtheta[3];
  real_T c2_upphi_theta;
  real_T c2_dx1_dpsi[3];
  real_T c2_dx2_dpsi[3];
  real_T c2_upphi_psi;
  real_T c2_delta1;
  real_T c2_delta2;
  real_T c2_delta3;
  real_T c2_ddot[3];
  real_T c2_nargin = 3.0;
  real_T c2_nargout = 6.0;
  real_T c2_x_ddot;
  real_T c2_y_ddot;
  real_T c2_z_ddot;
  real_T c2_phi_ddot;
  real_T c2_theta_ddot;
  real_T c2_psi_ddot;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_a;
  real_T c2_b;
  real_T c2_b_y;
  real_T c2_d_x;
  real_T c2_e_x;
  real_T c2_b_a;
  real_T c2_b_b;
  real_T c2_c_y;
  real_T c2_f_x;
  real_T c2_g_x;
  real_T c2_c_a;
  real_T c2_c_b;
  real_T c2_d_y;
  real_T c2_h_x;
  real_T c2_i_x;
  real_T c2_d_a;
  real_T c2_d_b;
  real_T c2_e_y;
  real_T c2_j_x;
  real_T c2_k_x;
  real_T c2_e_a;
  real_T c2_e_b;
  real_T c2_f_y;
  real_T c2_l_x;
  real_T c2_m_x;
  real_T c2_f_a;
  real_T c2_f_b;
  real_T c2_g_y;
  real_T c2_n_x;
  real_T c2_o_x;
  real_T c2_g_a;
  real_T c2_g_b;
  real_T c2_h_y;
  real_T c2_p_x;
  real_T c2_q_x;
  real_T c2_h_a;
  real_T c2_h_b;
  real_T c2_i_y;
  real_T c2_r_x;
  real_T c2_s_x;
  real_T c2_i_a;
  real_T c2_i_b;
  real_T c2_j_y;
  real_T c2_t_x;
  real_T c2_u_x;
  real_T c2_j_a;
  real_T c2_j_b;
  real_T c2_k_y;
  real_T c2_v_x;
  real_T c2_w_x;
  real_T c2_k_a;
  real_T c2_k_b;
  real_T c2_l_y;
  real_T c2_x_x;
  real_T c2_y_x;
  real_T c2_l_a;
  real_T c2_l_b;
  real_T c2_m_y;
  real_T c2_ab_x;
  real_T c2_bb_x;
  real_T c2_m_a;
  real_T c2_m_b;
  real_T c2_n_y;
  real_T c2_cb_x;
  real_T c2_db_x;
  real_T c2_n_a;
  real_T c2_n_b;
  real_T c2_o_y;
  real_T c2_eb_x;
  real_T c2_fb_x;
  real_T c2_o_a;
  real_T c2_o_b;
  real_T c2_p_y;
  real_T c2_gb_x;
  real_T c2_hb_x;
  real_T c2_p_a;
  real_T c2_p_b;
  real_T c2_q_y;
  real_T c2_ib_x;
  real_T c2_jb_x;
  real_T c2_q_a;
  real_T c2_q_b;
  real_T c2_r_y;
  real_T c2_kb_x;
  real_T c2_lb_x;
  real_T c2_r_a;
  real_T c2_r_b;
  real_T c2_s_y;
  real_T c2_mb_x;
  real_T c2_nb_x;
  real_T c2_s_a;
  real_T c2_s_b;
  real_T c2_t_y;
  real_T c2_ob_x;
  real_T c2_pb_x;
  real_T c2_t_a;
  real_T c2_t_b;
  real_T c2_u_y;
  real_T c2_qb_x;
  real_T c2_rb_x;
  real_T c2_u_a;
  real_T c2_u_b;
  real_T c2_v_y;
  real_T c2_v_a;
  real_T c2_v_b;
  real_T c2_w_y;
  real_T c2_w_a;
  real_T c2_w_b;
  real_T c2_x_y;
  real_T c2_x_a;
  real_T c2_x_b;
  real_T c2_y_y;
  real_T c2_y_a;
  real_T c2_y_b;
  real_T c2_ab_y;
  real_T c2_sb_x;
  real_T c2_tb_x;
  real_T c2_ab_a;
  real_T c2_ab_b;
  real_T c2_bb_y;
  real_T c2_ub_x;
  real_T c2_vb_x;
  real_T c2_bb_a;
  real_T c2_bb_b;
  real_T c2_cb_y;
  real_T c2_cb_a;
  real_T c2_cb_b;
  real_T c2_db_y;
  real_T c2_db_a;
  real_T c2_db_b;
  real_T c2_eb_y;
  real_T c2_eb_a;
  real_T c2_eb_b;
  real_T c2_fb_y;
  real_T c2_fb_a;
  real_T c2_fb_b;
  real_T c2_gb_y;
  real_T c2_wb_x;
  real_T c2_xb_x;
  real_T c2_yb_x;
  real_T c2_ac_x;
  real_T c2_gb_a;
  real_T c2_gb_b;
  real_T c2_hb_y;
  real_T c2_bc_x;
  real_T c2_cc_x;
  real_T c2_dc_x;
  real_T c2_ec_x;
  real_T c2_hb_a;
  real_T c2_hb_b;
  real_T c2_ib_y;
  real_T c2_fc_x;
  real_T c2_gc_x;
  real_T c2_ib_a;
  real_T c2_ib_b;
  real_T c2_jb_y;
  real_T c2_jb_a;
  real_T c2_jb_b;
  real_T c2_kb_y;
  real_T c2_kb_a;
  real_T c2_kb_b;
  real_T c2_lb_y;
  real_T c2_lb_a;
  real_T c2_lb_b;
  real_T c2_mb_y;
  real_T c2_mb_a;
  real_T c2_mb_b;
  real_T c2_nb_y;
  real_T c2_nb_a;
  real_T c2_nb_b;
  real_T c2_ob_y;
  real_T c2_hc_x;
  real_T c2_ic_x;
  real_T c2_jc_x;
  real_T c2_kc_x;
  real_T c2_ob_a;
  real_T c2_ob_b;
  real_T c2_pb_y;
  real_T c2_lc_x;
  real_T c2_mc_x;
  real_T c2_nc_x;
  real_T c2_oc_x;
  real_T c2_pb_a;
  real_T c2_pb_b;
  real_T c2_qb_y;
  real_T c2_pc_x;
  real_T c2_qc_x;
  real_T c2_qb_a;
  real_T c2_qb_b;
  real_T c2_rb_y;
  real_T c2_rb_a;
  real_T c2_rb_b;
  real_T c2_sb_y;
  real_T c2_sb_a;
  real_T c2_sb_b;
  real_T c2_tb_y;
  real_T c2_tb_a;
  real_T c2_tb_b;
  real_T c2_ub_y;
  real_T c2_ub_a;
  real_T c2_ub_b;
  real_T c2_vb_y;
  real_T c2_vb_a;
  real_T c2_vb_b;
  real_T c2_wb_y;
  real_T c2_rc_x;
  real_T c2_sc_x;
  real_T c2_wb_a;
  real_T c2_wb_b;
  real_T c2_xb_y;
  real_T c2_tc_x;
  real_T c2_uc_x;
  real_T c2_xb_a;
  real_T c2_xb_b;
  real_T c2_yb_y;
  real_T c2_yb_a;
  real_T c2_yb_b;
  real_T c2_ac_y;
  real_T c2_ac_a;
  real_T c2_ac_b;
  real_T c2_bc_y;
  real_T c2_bc_a;
  real_T c2_bc_b;
  real_T c2_cc_y;
  real_T c2_cc_a;
  real_T c2_cc_b;
  real_T c2_dc_y;
  real_T c2_vc_x;
  real_T c2_wc_x;
  real_T c2_xc_x;
  real_T c2_yc_x;
  real_T c2_dc_a;
  real_T c2_dc_b;
  real_T c2_ec_y;
  real_T c2_ad_x;
  real_T c2_bd_x;
  real_T c2_cd_x;
  real_T c2_dd_x;
  real_T c2_ec_a;
  real_T c2_ec_b;
  real_T c2_fc_y;
  real_T c2_ed_x;
  real_T c2_fd_x;
  real_T c2_fc_a;
  real_T c2_fc_b;
  real_T c2_gc_y;
  real_T c2_gc_a;
  real_T c2_gc_b;
  real_T c2_hc_y;
  real_T c2_hc_a;
  real_T c2_hc_b;
  real_T c2_ic_y;
  real_T c2_ic_a;
  real_T c2_ic_b;
  real_T c2_jc_y;
  real_T c2_jc_a;
  real_T c2_jc_b;
  real_T c2_kc_y;
  real_T c2_kc_a;
  real_T c2_kc_b;
  real_T c2_lc_y;
  real_T c2_gd_x;
  real_T c2_hd_x;
  real_T c2_id_x;
  real_T c2_jd_x;
  real_T c2_lc_a;
  real_T c2_lc_b;
  real_T c2_mc_y;
  real_T c2_kd_x;
  real_T c2_ld_x;
  real_T c2_md_x;
  real_T c2_nd_x;
  real_T c2_mc_a;
  real_T c2_mc_b;
  real_T c2_nc_y;
  real_T c2_od_x;
  real_T c2_pd_x;
  real_T c2_nc_a;
  real_T c2_nc_b;
  real_T c2_oc_y;
  real_T c2_oc_a;
  real_T c2_oc_b;
  real_T c2_pc_y;
  real_T c2_pc_a;
  real_T c2_pc_b;
  real_T c2_qc_y;
  real_T c2_qc_a;
  real_T c2_qc_b;
  real_T c2_rc_y;
  real_T c2_rc_a;
  real_T c2_rc_b;
  real_T c2_sc_y;
  real_T c2_sc_a;
  real_T c2_sc_b;
  real_T c2_tc_y;
  real_T c2_qd_x;
  real_T c2_rd_x;
  real_T c2_tc_a;
  real_T c2_tc_b;
  real_T c2_uc_y;
  real_T c2_uc_a;
  real_T c2_uc_b;
  real_T c2_vc_y;
  real_T c2_vc_a;
  real_T c2_vc_b;
  real_T c2_wc_y;
  real_T c2_wc_a;
  real_T c2_wc_b;
  real_T c2_xc_y;
  real_T c2_xc_a;
  real_T c2_xc_b;
  real_T c2_yc_y;
  real_T c2_sd_x;
  real_T c2_td_x;
  real_T c2_yc_a;
  real_T c2_yc_b;
  real_T c2_ad_y;
  real_T c2_ud_x;
  real_T c2_vd_x;
  real_T c2_ad_a;
  real_T c2_ad_b;
  real_T c2_bd_y;
  real_T c2_bd_a;
  real_T c2_bd_b;
  real_T c2_cd_y;
  real_T c2_cd_a;
  real_T c2_cd_b;
  real_T c2_dd_y;
  real_T c2_dd_a;
  real_T c2_dd_b;
  real_T c2_ed_y;
  real_T c2_ed_a;
  real_T c2_ed_b;
  real_T c2_fd_y;
  real_T c2_wd_x;
  real_T c2_xd_x;
  real_T c2_fd_a;
  real_T c2_fd_b;
  real_T c2_gd_y;
  real_T c2_yd_x;
  real_T c2_ae_x;
  real_T c2_gd_a;
  real_T c2_gd_b;
  real_T c2_hd_y;
  real_T c2_be_x;
  real_T c2_ce_x;
  real_T c2_de_x;
  real_T c2_ee_x;
  real_T c2_hd_a;
  real_T c2_hd_b;
  real_T c2_id_y;
  real_T c2_fe_x;
  real_T c2_ge_x;
  real_T c2_he_x;
  real_T c2_ie_x;
  real_T c2_id_a;
  real_T c2_id_b;
  real_T c2_jd_y;
  real_T c2_je_x;
  real_T c2_ke_x;
  real_T c2_jd_a;
  real_T c2_jd_b;
  real_T c2_kd_y;
  real_T c2_le_x;
  real_T c2_me_x;
  real_T c2_ne_x;
  real_T c2_oe_x;
  real_T c2_kd_a;
  real_T c2_kd_b;
  real_T c2_ld_y;
  real_T c2_pe_x;
  real_T c2_qe_x;
  real_T c2_re_x;
  real_T c2_se_x;
  real_T c2_ld_a;
  real_T c2_ld_b;
  real_T c2_md_y;
  real_T c2_te_x;
  real_T c2_ue_x;
  real_T c2_md_a;
  real_T c2_md_b;
  real_T c2_nd_y;
  real_T c2_ve_x;
  real_T c2_we_x;
  real_T c2_xe_x;
  real_T c2_ye_x;
  real_T c2_nd_a;
  real_T c2_nd_b;
  real_T c2_od_y;
  real_T c2_af_x;
  real_T c2_bf_x;
  real_T c2_cf_x;
  real_T c2_df_x;
  real_T c2_od_a;
  real_T c2_od_b;
  real_T c2_pd_y;
  real_T c2_ef_x;
  real_T c2_ff_x;
  real_T c2_gf_x;
  real_T c2_hf_x;
  real_T c2_pd_a;
  real_T c2_pd_b;
  real_T c2_qd_y;
  real_T c2_if_x;
  real_T c2_jf_x;
  real_T c2_qd_a;
  real_T c2_qd_b;
  real_T c2_rd_y;
  real_T c2_kf_x;
  real_T c2_lf_x;
  real_T c2_mf_x;
  real_T c2_nf_x;
  real_T c2_rd_a;
  real_T c2_rd_b;
  real_T c2_sd_y;
  real_T c2_of_x;
  real_T c2_pf_x;
  real_T c2_qf_x;
  real_T c2_rf_x;
  real_T c2_sd_a;
  real_T c2_sd_b;
  real_T c2_td_y;
  real_T c2_sf_x;
  real_T c2_tf_x;
  real_T c2_td_a;
  real_T c2_td_b;
  real_T c2_ud_y;
  real_T c2_uf_x;
  real_T c2_vf_x;
  real_T c2_wf_x;
  real_T c2_xf_x;
  real_T c2_ud_a;
  real_T c2_ud_b;
  real_T c2_vd_y;
  real_T c2_yf_x;
  real_T c2_ag_x;
  real_T c2_bg_x;
  real_T c2_cg_x;
  real_T c2_vd_a;
  real_T c2_vd_b;
  real_T c2_wd_y;
  real_T c2_dg_x;
  real_T c2_eg_x;
  real_T c2_fg_x;
  real_T c2_gg_x;
  real_T c2_wd_a;
  real_T c2_wd_b;
  real_T c2_xd_y;
  real_T c2_hg_x;
  real_T c2_ig_x;
  real_T c2_xd_a;
  real_T c2_xd_b;
  real_T c2_yd_y;
  real_T c2_yd_a;
  real_T c2_yd_b;
  real_T c2_ae_y;
  real_T c2_ae_a;
  real_T c2_ae_b;
  real_T c2_be_y;
  real_T c2_be_a;
  real_T c2_be_b;
  real_T c2_ce_y;
  real_T c2_ce_a;
  real_T c2_ce_b;
  real_T c2_de_y;
  real_T c2_de_a;
  real_T c2_de_b;
  real_T c2_ee_y;
  real_T c2_ee_a;
  real_T c2_ee_b;
  real_T c2_fe_y;
  real_T c2_fe_a;
  real_T c2_fe_b;
  real_T c2_ge_y;
  real_T c2_ge_a;
  real_T c2_ge_b;
  real_T c2_he_y;
  real_T c2_he_a;
  real_T c2_he_b;
  real_T c2_ie_y;
  real_T c2_ie_a;
  real_T c2_ie_b;
  real_T c2_je_y;
  real_T c2_je_a;
  real_T c2_je_b;
  real_T c2_ke_y;
  real_T c2_A;
  real_T c2_jg_x;
  real_T c2_kg_x;
  real_T c2_le_y;
  real_T c2_lg_x;
  real_T c2_mg_x;
  real_T c2_ng_x;
  real_T c2_og_x;
  real_T c2_ke_a;
  real_T c2_ke_b;
  real_T c2_me_y;
  real_T c2_pg_x;
  real_T c2_qg_x;
  real_T c2_rg_x;
  real_T c2_sg_x;
  real_T c2_le_a;
  real_T c2_le_b;
  real_T c2_ne_y;
  real_T c2_tg_x;
  real_T c2_ug_x;
  real_T c2_me_a;
  real_T c2_me_b;
  real_T c2_oe_y;
  real_T c2_vg_x;
  real_T c2_wg_x;
  real_T c2_xg_x;
  real_T c2_yg_x;
  real_T c2_ne_a;
  real_T c2_ne_b;
  real_T c2_pe_y;
  real_T c2_ah_x;
  real_T c2_bh_x;
  real_T c2_ch_x;
  real_T c2_dh_x;
  real_T c2_oe_a;
  real_T c2_oe_b;
  real_T c2_qe_y;
  real_T c2_eh_x;
  real_T c2_fh_x;
  real_T c2_pe_a;
  real_T c2_pe_b;
  real_T c2_re_y;
  real_T c2_gh_x;
  real_T c2_hh_x;
  real_T c2_ih_x;
  real_T c2_jh_x;
  real_T c2_qe_a;
  real_T c2_qe_b;
  real_T c2_se_y;
  real_T c2_kh_x;
  real_T c2_lh_x;
  real_T c2_mh_x;
  real_T c2_nh_x;
  real_T c2_re_a;
  real_T c2_re_b;
  real_T c2_te_y;
  real_T c2_oh_x;
  real_T c2_ph_x;
  real_T c2_se_a;
  real_T c2_se_b;
  real_T c2_ue_y;
  real_T c2_qh_x;
  real_T c2_rh_x;
  real_T c2_sh_x;
  real_T c2_th_x;
  real_T c2_te_a;
  real_T c2_te_b;
  real_T c2_ve_y;
  real_T c2_uh_x;
  real_T c2_vh_x;
  real_T c2_wh_x;
  real_T c2_xh_x;
  real_T c2_ue_a;
  real_T c2_ue_b;
  real_T c2_we_y;
  real_T c2_yh_x;
  real_T c2_ai_x;
  real_T c2_ve_a;
  real_T c2_ve_b;
  real_T c2_xe_y;
  real_T c2_we_a;
  real_T c2_we_b[3];
  int32_T c2_i6;
  real_T c2_b_A;
  real_T c2_bi_x;
  real_T c2_ci_x;
  real_T c2_ye_y;
  real_T c2_di_x;
  real_T c2_ei_x;
  real_T c2_fi_x;
  real_T c2_gi_x;
  real_T c2_xe_a;
  real_T c2_xe_b;
  real_T c2_af_y;
  real_T c2_hi_x;
  real_T c2_ii_x;
  real_T c2_ji_x;
  real_T c2_ki_x;
  real_T c2_ye_a;
  real_T c2_ye_b;
  real_T c2_bf_y;
  real_T c2_li_x;
  real_T c2_mi_x;
  real_T c2_af_a;
  real_T c2_af_b;
  real_T c2_cf_y;
  real_T c2_ni_x;
  real_T c2_oi_x;
  real_T c2_pi_x;
  real_T c2_qi_x;
  real_T c2_bf_a;
  real_T c2_bf_b;
  real_T c2_df_y;
  real_T c2_ri_x;
  real_T c2_si_x;
  real_T c2_ti_x;
  real_T c2_ui_x;
  real_T c2_cf_a;
  real_T c2_cf_b;
  real_T c2_ef_y;
  real_T c2_vi_x;
  real_T c2_wi_x;
  real_T c2_df_a;
  real_T c2_df_b;
  real_T c2_ff_y;
  real_T c2_xi_x;
  real_T c2_yi_x;
  real_T c2_aj_x;
  real_T c2_bj_x;
  real_T c2_ef_a;
  real_T c2_ef_b;
  real_T c2_gf_y;
  real_T c2_cj_x;
  real_T c2_dj_x;
  real_T c2_ej_x;
  real_T c2_fj_x;
  real_T c2_ff_a;
  real_T c2_ff_b;
  real_T c2_hf_y;
  real_T c2_gj_x;
  real_T c2_hj_x;
  real_T c2_gf_a;
  real_T c2_gf_b;
  real_T c2_if_y;
  real_T c2_ij_x;
  real_T c2_jj_x;
  real_T c2_kj_x;
  real_T c2_lj_x;
  real_T c2_hf_a;
  real_T c2_hf_b;
  real_T c2_jf_y;
  real_T c2_mj_x;
  real_T c2_nj_x;
  real_T c2_oj_x;
  real_T c2_pj_x;
  real_T c2_if_a;
  real_T c2_if_b;
  real_T c2_kf_y;
  real_T c2_qj_x;
  real_T c2_rj_x;
  real_T c2_jf_a;
  real_T c2_jf_b;
  real_T c2_lf_y;
  real_T c2_kf_a;
  int32_T c2_i7;
  int32_T c2_i8;
  real_T c2_lf_a[9];
  int32_T c2_i9;
  int32_T c2_i10;
  real_T c2_mf_y[3];
  int32_T c2_i11;
  int32_T c2_i12;
  int32_T c2_i13;
  int32_T c2_i14;
  int32_T c2_i15;
  real_T c2_nf_y[3];
  int32_T c2_i16;
  int32_T c2_i17;
  int32_T c2_i18;
  real_T c2_of_y[3];
  int32_T c2_i19;
  real_T c2_b_dx1_dphi[3];
  int32_T c2_i20;
  real_T c2_pf_y[3];
  int32_T c2_i21;
  real_T c2_b_dx2_dphi[3];
  real_T c2_c_A;
  real_T c2_sj_x;
  real_T c2_tj_x;
  real_T c2_qf_y;
  real_T c2_uj_x;
  real_T c2_vj_x;
  real_T c2_wj_x;
  real_T c2_xj_x;
  real_T c2_mf_a;
  real_T c2_kf_b;
  real_T c2_rf_y;
  real_T c2_yj_x;
  real_T c2_ak_x;
  real_T c2_nf_a;
  real_T c2_lf_b;
  real_T c2_sf_y;
  real_T c2_bk_x;
  real_T c2_ck_x;
  real_T c2_dk_x;
  real_T c2_ek_x;
  real_T c2_of_a;
  real_T c2_mf_b;
  real_T c2_tf_y;
  real_T c2_fk_x;
  real_T c2_gk_x;
  real_T c2_pf_a;
  real_T c2_nf_b;
  real_T c2_uf_y;
  real_T c2_hk_x;
  real_T c2_ik_x;
  real_T c2_jk_x;
  real_T c2_kk_x;
  real_T c2_qf_a;
  real_T c2_of_b;
  real_T c2_vf_y;
  real_T c2_lk_x;
  real_T c2_mk_x;
  real_T c2_rf_a;
  real_T c2_pf_b;
  real_T c2_wf_y;
  real_T c2_nk_x;
  real_T c2_ok_x;
  real_T c2_pk_x;
  real_T c2_qk_x;
  real_T c2_sf_a;
  real_T c2_qf_b;
  real_T c2_xf_y;
  real_T c2_rk_x;
  real_T c2_sk_x;
  real_T c2_tf_a;
  real_T c2_rf_b;
  real_T c2_yf_y;
  real_T c2_tk_x;
  real_T c2_uk_x;
  real_T c2_vk_x;
  real_T c2_wk_x;
  real_T c2_uf_a;
  real_T c2_sf_b;
  real_T c2_ag_y;
  real_T c2_xk_x;
  real_T c2_yk_x;
  real_T c2_al_x;
  real_T c2_bl_x;
  real_T c2_vf_a;
  real_T c2_tf_b;
  real_T c2_bg_y;
  real_T c2_wf_a;
  int32_T c2_i22;
  real_T c2_d_A;
  real_T c2_cl_x;
  real_T c2_dl_x;
  real_T c2_cg_y;
  real_T c2_el_x;
  real_T c2_fl_x;
  real_T c2_gl_x;
  real_T c2_hl_x;
  real_T c2_xf_a;
  real_T c2_uf_b;
  real_T c2_dg_y;
  real_T c2_il_x;
  real_T c2_jl_x;
  real_T c2_yf_a;
  real_T c2_vf_b;
  real_T c2_eg_y;
  real_T c2_kl_x;
  real_T c2_ll_x;
  real_T c2_ml_x;
  real_T c2_nl_x;
  real_T c2_ag_a;
  real_T c2_wf_b;
  real_T c2_fg_y;
  real_T c2_ol_x;
  real_T c2_pl_x;
  real_T c2_bg_a;
  real_T c2_xf_b;
  real_T c2_gg_y;
  real_T c2_ql_x;
  real_T c2_rl_x;
  real_T c2_sl_x;
  real_T c2_tl_x;
  real_T c2_cg_a;
  real_T c2_yf_b;
  real_T c2_hg_y;
  real_T c2_ul_x;
  real_T c2_vl_x;
  real_T c2_dg_a;
  real_T c2_ag_b;
  real_T c2_ig_y;
  real_T c2_wl_x;
  real_T c2_xl_x;
  real_T c2_yl_x;
  real_T c2_am_x;
  real_T c2_eg_a;
  real_T c2_bg_b;
  real_T c2_jg_y;
  real_T c2_bm_x;
  real_T c2_cm_x;
  real_T c2_fg_a;
  real_T c2_cg_b;
  real_T c2_kg_y;
  real_T c2_dm_x;
  real_T c2_em_x;
  real_T c2_fm_x;
  real_T c2_gm_x;
  real_T c2_gg_a;
  real_T c2_dg_b;
  real_T c2_lg_y;
  real_T c2_hm_x;
  real_T c2_im_x;
  real_T c2_jm_x;
  real_T c2_km_x;
  real_T c2_hg_a;
  real_T c2_eg_b;
  real_T c2_mg_y;
  real_T c2_ig_a;
  int32_T c2_i23;
  int32_T c2_i24;
  int32_T c2_i25;
  int32_T c2_i26;
  int32_T c2_i27;
  int32_T c2_i28;
  int32_T c2_i29;
  int32_T c2_i30;
  int32_T c2_i31;
  int32_T c2_i32;
  int32_T c2_i33;
  int32_T c2_i34;
  real_T c2_ng_y[3];
  int32_T c2_i35;
  real_T c2_b_dx1_dtheta[3];
  int32_T c2_i36;
  real_T c2_og_y[3];
  int32_T c2_i37;
  real_T c2_b_dx2_dtheta[3];
  real_T c2_e_A;
  real_T c2_lm_x;
  real_T c2_mm_x;
  real_T c2_pg_y;
  real_T c2_nm_x;
  real_T c2_om_x;
  real_T c2_pm_x;
  real_T c2_qm_x;
  real_T c2_jg_a;
  real_T c2_fg_b;
  real_T c2_qg_y;
  real_T c2_rm_x;
  real_T c2_sm_x;
  real_T c2_tm_x;
  real_T c2_um_x;
  real_T c2_kg_a;
  real_T c2_gg_b;
  real_T c2_rg_y;
  real_T c2_vm_x;
  real_T c2_wm_x;
  real_T c2_lg_a;
  real_T c2_hg_b;
  real_T c2_sg_y;
  real_T c2_xm_x;
  real_T c2_ym_x;
  real_T c2_an_x;
  real_T c2_bn_x;
  real_T c2_mg_a;
  real_T c2_ig_b;
  real_T c2_tg_y;
  real_T c2_cn_x;
  real_T c2_dn_x;
  real_T c2_en_x;
  real_T c2_fn_x;
  real_T c2_ng_a;
  real_T c2_jg_b;
  real_T c2_ug_y;
  real_T c2_gn_x;
  real_T c2_hn_x;
  real_T c2_og_a;
  real_T c2_kg_b;
  real_T c2_vg_y;
  real_T c2_in_x;
  real_T c2_jn_x;
  real_T c2_kn_x;
  real_T c2_ln_x;
  real_T c2_pg_a;
  real_T c2_lg_b;
  real_T c2_wg_y;
  real_T c2_mn_x;
  real_T c2_nn_x;
  real_T c2_on_x;
  real_T c2_pn_x;
  real_T c2_qg_a;
  real_T c2_mg_b;
  real_T c2_xg_y;
  real_T c2_qn_x;
  real_T c2_rn_x;
  real_T c2_rg_a;
  real_T c2_ng_b;
  real_T c2_yg_y;
  real_T c2_sn_x;
  real_T c2_tn_x;
  real_T c2_un_x;
  real_T c2_vn_x;
  real_T c2_sg_a;
  real_T c2_og_b;
  real_T c2_ah_y;
  real_T c2_wn_x;
  real_T c2_xn_x;
  real_T c2_yn_x;
  real_T c2_ao_x;
  real_T c2_tg_a;
  real_T c2_pg_b;
  real_T c2_bh_y;
  real_T c2_bo_x;
  real_T c2_co_x;
  real_T c2_ug_a;
  real_T c2_qg_b;
  real_T c2_ch_y;
  real_T c2_do_x;
  real_T c2_eo_x;
  real_T c2_fo_x;
  real_T c2_go_x;
  real_T c2_vg_a;
  real_T c2_rg_b;
  real_T c2_dh_y;
  real_T c2_ho_x;
  real_T c2_io_x;
  real_T c2_jo_x;
  real_T c2_ko_x;
  real_T c2_wg_a;
  real_T c2_sg_b;
  real_T c2_eh_y;
  real_T c2_xg_a;
  int32_T c2_i38;
  real_T c2_f_A;
  real_T c2_lo_x;
  real_T c2_mo_x;
  real_T c2_fh_y;
  real_T c2_no_x;
  real_T c2_oo_x;
  real_T c2_po_x;
  real_T c2_qo_x;
  real_T c2_yg_a;
  real_T c2_tg_b;
  real_T c2_gh_y;
  real_T c2_ro_x;
  real_T c2_so_x;
  real_T c2_to_x;
  real_T c2_uo_x;
  real_T c2_ah_a;
  real_T c2_ug_b;
  real_T c2_hh_y;
  real_T c2_vo_x;
  real_T c2_wo_x;
  real_T c2_bh_a;
  real_T c2_vg_b;
  real_T c2_ih_y;
  real_T c2_xo_x;
  real_T c2_yo_x;
  real_T c2_ap_x;
  real_T c2_bp_x;
  real_T c2_ch_a;
  real_T c2_wg_b;
  real_T c2_jh_y;
  real_T c2_cp_x;
  real_T c2_dp_x;
  real_T c2_ep_x;
  real_T c2_fp_x;
  real_T c2_dh_a;
  real_T c2_xg_b;
  real_T c2_kh_y;
  real_T c2_gp_x;
  real_T c2_hp_x;
  real_T c2_eh_a;
  real_T c2_yg_b;
  real_T c2_lh_y;
  real_T c2_ip_x;
  real_T c2_jp_x;
  real_T c2_kp_x;
  real_T c2_lp_x;
  real_T c2_fh_a;
  real_T c2_ah_b;
  real_T c2_mh_y;
  real_T c2_mp_x;
  real_T c2_np_x;
  real_T c2_op_x;
  real_T c2_pp_x;
  real_T c2_gh_a;
  real_T c2_bh_b;
  real_T c2_nh_y;
  real_T c2_qp_x;
  real_T c2_rp_x;
  real_T c2_hh_a;
  real_T c2_ch_b;
  real_T c2_oh_y;
  real_T c2_sp_x;
  real_T c2_tp_x;
  real_T c2_up_x;
  real_T c2_vp_x;
  real_T c2_ih_a;
  real_T c2_dh_b;
  real_T c2_ph_y;
  real_T c2_wp_x;
  real_T c2_xp_x;
  real_T c2_yp_x;
  real_T c2_aq_x;
  real_T c2_jh_a;
  real_T c2_eh_b;
  real_T c2_qh_y;
  real_T c2_bq_x;
  real_T c2_cq_x;
  real_T c2_kh_a;
  real_T c2_fh_b;
  real_T c2_rh_y;
  real_T c2_dq_x;
  real_T c2_eq_x;
  real_T c2_fq_x;
  real_T c2_gq_x;
  real_T c2_lh_a;
  real_T c2_gh_b;
  real_T c2_sh_y;
  real_T c2_hq_x;
  real_T c2_iq_x;
  real_T c2_jq_x;
  real_T c2_kq_x;
  real_T c2_mh_a;
  real_T c2_hh_b;
  real_T c2_th_y;
  real_T c2_nh_a;
  int32_T c2_i39;
  int32_T c2_i40;
  int32_T c2_i41;
  int32_T c2_i42;
  int32_T c2_i43;
  int32_T c2_i44;
  int32_T c2_i45;
  int32_T c2_i46;
  int32_T c2_i47;
  int32_T c2_i48;
  int32_T c2_i49;
  int32_T c2_i50;
  real_T c2_uh_y[3];
  int32_T c2_i51;
  real_T c2_b_dx1_dpsi[3];
  int32_T c2_i52;
  real_T c2_vh_y[3];
  int32_T c2_i53;
  real_T c2_b_dx2_dpsi[3];
  real_T c2_ih_b;
  real_T c2_wh_y;
  real_T c2_lq_x;
  real_T c2_mq_x;
  real_T c2_oh_a;
  real_T c2_jh_b;
  real_T c2_xh_y;
  real_T c2_nq_x;
  real_T c2_oq_x;
  real_T c2_ph_a;
  real_T c2_kh_b;
  real_T c2_yh_y;
  real_T c2_pq_x;
  real_T c2_qq_x;
  real_T c2_qh_a;
  real_T c2_lh_b;
  real_T c2_ai_y;
  real_T c2_rh_a;
  real_T c2_mh_b;
  real_T c2_bi_y;
  real_T c2_sh_a;
  real_T c2_nh_b;
  real_T c2_ci_y;
  real_T c2_oh_b;
  real_T c2_di_y;
  real_T c2_rq_x;
  real_T c2_sq_x;
  real_T c2_th_a;
  real_T c2_ph_b;
  real_T c2_ei_y;
  real_T c2_tq_x;
  real_T c2_uq_x;
  real_T c2_uh_a;
  real_T c2_qh_b;
  real_T c2_fi_y;
  real_T c2_vq_x;
  real_T c2_wq_x;
  real_T c2_vh_a;
  real_T c2_rh_b;
  real_T c2_gi_y;
  real_T c2_wh_a;
  real_T c2_sh_b;
  real_T c2_hi_y;
  real_T c2_xh_a;
  real_T c2_th_b;
  real_T c2_ii_y;
  real_T c2_xq_x;
  real_T c2_yq_x;
  real_T c2_yh_a;
  real_T c2_uh_b;
  real_T c2_ji_y;
  real_T c2_ar_x;
  real_T c2_br_x;
  real_T c2_ai_a;
  real_T c2_vh_b;
  real_T c2_ki_y;
  real_T c2_bi_a;
  real_T c2_wh_b;
  real_T c2_li_y;
  real_T c2_ci_a;
  real_T c2_xh_b;
  real_T c2_mi_y;
  real_T c2_cr_x;
  real_T c2_dr_x;
  real_T c2_di_a;
  real_T c2_yh_b;
  real_T c2_ni_y;
  real_T c2_er_x;
  real_T c2_fr_x;
  real_T c2_ei_a;
  real_T c2_ai_b;
  real_T c2_oi_y;
  real_T c2_gr_x;
  real_T c2_hr_x;
  real_T c2_fi_a;
  real_T c2_bi_b;
  real_T c2_pi_y;
  real_T c2_gi_a;
  real_T c2_ci_b;
  real_T c2_qi_y;
  real_T c2_ir_x;
  real_T c2_jr_x;
  real_T c2_hi_a;
  real_T c2_di_b;
  real_T c2_ri_y;
  real_T c2_kr_x;
  real_T c2_lr_x;
  real_T c2_ii_a;
  real_T c2_ei_b;
  real_T c2_si_y;
  real_T c2_ji_a;
  real_T c2_fi_b;
  real_T c2_ti_y;
  real_T c2_ki_a;
  real_T c2_gi_b;
  real_T c2_ui_y;
  real_T c2_hi_b;
  real_T c2_vi_y;
  real_T c2_mr_x;
  real_T c2_nr_x;
  real_T c2_li_a;
  real_T c2_ii_b;
  real_T c2_wi_y;
  real_T c2_or_x;
  real_T c2_pr_x;
  real_T c2_mi_a;
  real_T c2_ji_b;
  real_T c2_xi_y;
  real_T c2_qr_x;
  real_T c2_rr_x;
  real_T c2_ni_a;
  real_T c2_ki_b;
  real_T c2_yi_y;
  real_T c2_oi_a;
  real_T c2_li_b;
  real_T c2_aj_y;
  real_T c2_pi_a;
  real_T c2_mi_b;
  real_T c2_bj_y;
  real_T c2_ni_b;
  real_T c2_cj_y;
  real_T c2_sr_x;
  real_T c2_tr_x;
  real_T c2_qi_a;
  real_T c2_oi_b;
  real_T c2_dj_y;
  real_T c2_ur_x;
  real_T c2_vr_x;
  real_T c2_ri_a;
  real_T c2_pi_b;
  real_T c2_ej_y;
  real_T c2_wr_x;
  real_T c2_xr_x;
  real_T c2_si_a;
  real_T c2_qi_b;
  real_T c2_fj_y;
  real_T c2_ti_a;
  real_T c2_ri_b;
  real_T c2_gj_y;
  real_T c2_ui_a;
  real_T c2_si_b;
  real_T c2_hj_y;
  real_T c2_yr_x;
  real_T c2_as_x;
  real_T c2_vi_a;
  real_T c2_ti_b;
  real_T c2_ij_y;
  real_T c2_bs_x;
  real_T c2_cs_x;
  real_T c2_wi_a;
  real_T c2_ui_b;
  real_T c2_jj_y;
  real_T c2_xi_a;
  real_T c2_vi_b;
  real_T c2_kj_y;
  real_T c2_yi_a;
  real_T c2_wi_b;
  real_T c2_lj_y;
  real_T c2_ds_x;
  real_T c2_es_x;
  real_T c2_aj_a;
  real_T c2_xi_b;
  real_T c2_mj_y;
  real_T c2_fs_x;
  real_T c2_gs_x;
  real_T c2_bj_a;
  real_T c2_yi_b;
  real_T c2_nj_y;
  real_T c2_hs_x;
  real_T c2_is_x;
  real_T c2_cj_a;
  real_T c2_aj_b;
  real_T c2_oj_y;
  real_T c2_dj_a;
  real_T c2_bj_b;
  real_T c2_pj_y;
  real_T c2_js_x;
  real_T c2_ks_x;
  real_T c2_ej_a;
  real_T c2_cj_b;
  real_T c2_qj_y;
  real_T c2_ls_x;
  real_T c2_ms_x;
  real_T c2_fj_a;
  real_T c2_dj_b;
  real_T c2_rj_y;
  real_T c2_gj_a;
  real_T c2_ej_b;
  real_T c2_sj_y;
  real_T c2_hj_a;
  real_T c2_fj_b;
  real_T c2_tj_y;
  real_T c2_gj_b;
  real_T c2_uj_y;
  real_T c2_ns_x;
  real_T c2_os_x;
  real_T c2_ij_a;
  real_T c2_hj_b;
  real_T c2_vj_y;
  real_T c2_ps_x;
  real_T c2_qs_x;
  real_T c2_jj_a;
  real_T c2_ij_b;
  real_T c2_wj_y;
  real_T c2_kj_a;
  real_T c2_jj_b;
  real_T c2_xj_y;
  real_T c2_lj_a;
  real_T c2_kj_b;
  real_T c2_yj_y;
  real_T c2_rs_x;
  real_T c2_ss_x;
  real_T c2_mj_a;
  real_T c2_lj_b;
  real_T c2_ak_y;
  real_T c2_nj_a;
  real_T c2_mj_b;
  real_T c2_bk_y;
  real_T c2_oj_a;
  real_T c2_nj_b;
  real_T c2_ck_y;
  real_T c2_ts_x;
  real_T c2_us_x;
  real_T c2_pj_a;
  real_T c2_oj_b;
  real_T c2_dk_y;
  real_T c2_vs_x;
  real_T c2_ws_x;
  real_T c2_qj_a;
  real_T c2_pj_b;
  real_T c2_ek_y;
  real_T c2_rj_a;
  real_T c2_qj_b;
  real_T c2_fk_y;
  real_T c2_sj_a;
  real_T c2_rj_b;
  real_T c2_gk_y;
  real_T c2_xs_x;
  real_T c2_ys_x;
  real_T c2_tj_a;
  real_T c2_sj_b;
  real_T c2_hk_y;
  real_T c2_at_x;
  real_T c2_bt_x;
  real_T c2_uj_a;
  real_T c2_tj_b;
  real_T c2_ik_y;
  real_T c2_vj_a;
  real_T c2_uj_b;
  real_T c2_jk_y;
  real_T c2_wj_a;
  real_T c2_vj_b;
  real_T c2_kk_y;
  real_T c2_wj_b;
  real_T c2_lk_y;
  real_T c2_ct_x;
  real_T c2_dt_x;
  real_T c2_xj_a;
  real_T c2_xj_b;
  real_T c2_mk_y;
  real_T c2_et_x;
  real_T c2_ft_x;
  real_T c2_yj_a;
  real_T c2_yj_b;
  real_T c2_nk_y;
  real_T c2_ak_a;
  real_T c2_ak_b;
  real_T c2_ok_y;
  real_T c2_bk_a;
  real_T c2_bk_b;
  real_T c2_pk_y;
  real_T c2_gt_x;
  real_T c2_ht_x;
  real_T c2_ck_a;
  real_T c2_ck_b;
  real_T c2_qk_y;
  real_T c2_it_x;
  real_T c2_jt_x;
  real_T c2_dk_a;
  real_T c2_dk_b;
  real_T c2_rk_y;
  real_T c2_kt_x;
  real_T c2_lt_x;
  real_T c2_ek_a;
  real_T c2_ek_b;
  real_T c2_sk_y;
  real_T c2_fk_a;
  real_T c2_fk_b;
  real_T c2_tk_y;
  real_T c2_mt_x;
  real_T c2_nt_x;
  real_T c2_gk_a;
  real_T c2_gk_b;
  real_T c2_uk_y;
  real_T c2_ot_x;
  real_T c2_pt_x;
  real_T c2_hk_a;
  real_T c2_hk_b;
  real_T c2_vk_y;
  real_T c2_ik_a;
  real_T c2_ik_b;
  real_T c2_wk_y;
  real_T c2_jk_a;
  real_T c2_jk_b;
  real_T c2_xk_y;
  real_T c2_qt_x;
  real_T c2_rt_x;
  real_T c2_kk_a;
  real_T c2_kk_b;
  real_T c2_yk_y;
  real_T c2_st_x;
  real_T c2_tt_x;
  real_T c2_lk_a;
  real_T c2_lk_b;
  real_T c2_al_y;
  real_T c2_mk_a;
  real_T c2_mk_b;
  real_T c2_bl_y;
  real_T c2_nk_a;
  real_T c2_nk_b;
  real_T c2_cl_y;
  real_T c2_ok_b;
  real_T c2_dl_y;
  real_T c2_ut_x;
  real_T c2_vt_x;
  real_T c2_ok_a;
  real_T c2_pk_b;
  real_T c2_el_y;
  real_T c2_wt_x;
  real_T c2_xt_x;
  real_T c2_pk_a;
  real_T c2_qk_b;
  real_T c2_fl_y;
  real_T c2_qk_a;
  real_T c2_rk_b;
  real_T c2_gl_y;
  real_T c2_rk_a;
  real_T c2_sk_b;
  real_T c2_hl_y;
  real_T c2_yt_x;
  real_T c2_au_x;
  real_T c2_sk_a;
  real_T c2_tk_b;
  real_T c2_il_y;
  real_T c2_bu_x;
  real_T c2_cu_x;
  real_T c2_tk_a;
  real_T c2_uk_b;
  real_T c2_jl_y;
  real_T c2_du_x;
  real_T c2_eu_x;
  real_T c2_uk_a;
  real_T c2_vk_b;
  real_T c2_kl_y;
  real_T c2_vk_a;
  real_T c2_wk_b;
  real_T c2_ll_y;
  real_T c2_fu_x;
  real_T c2_gu_x;
  real_T c2_wk_a;
  real_T c2_xk_b;
  real_T c2_ml_y;
  real_T c2_hu_x;
  real_T c2_iu_x;
  real_T c2_xk_a;
  real_T c2_yk_b;
  real_T c2_nl_y;
  real_T c2_yk_a;
  real_T c2_al_b;
  real_T c2_ol_y;
  real_T c2_ju_x;
  real_T c2_ku_x;
  real_T c2_al_a;
  real_T c2_bl_b;
  real_T c2_pl_y;
  real_T c2_bl_a;
  real_T c2_cl_b;
  real_T c2_ql_y;
  real_T c2_cl_a;
  real_T c2_dl_b;
  real_T c2_rl_y;
  real_T c2_lu_x;
  real_T c2_mu_x;
  real_T c2_dl_a;
  real_T c2_el_b;
  real_T c2_sl_y;
  real_T c2_nu_x;
  real_T c2_ou_x;
  real_T c2_el_a;
  real_T c2_fl_b;
  real_T c2_tl_y;
  real_T c2_pu_x;
  real_T c2_qu_x;
  real_T c2_fl_a;
  real_T c2_gl_b;
  real_T c2_ul_y;
  real_T c2_gl_a;
  real_T c2_hl_b;
  real_T c2_vl_y;
  real_T c2_ru_x;
  real_T c2_su_x;
  real_T c2_hl_a;
  real_T c2_il_b;
  real_T c2_wl_y;
  real_T c2_tu_x;
  real_T c2_uu_x;
  real_T c2_il_a;
  real_T c2_jl_b;
  real_T c2_xl_y;
  real_T c2_jl_a;
  real_T c2_kl_b;
  real_T c2_yl_y;
  real_T c2_kl_a;
  real_T c2_ll_b;
  real_T c2_am_y;
  real_T c2_vu_x;
  real_T c2_wu_x;
  real_T c2_ll_a;
  real_T c2_ml_b;
  real_T c2_bm_y;
  real_T c2_xu_x;
  real_T c2_yu_x;
  real_T c2_ml_a;
  real_T c2_nl_b;
  real_T c2_cm_y;
  real_T c2_nl_a;
  real_T c2_ol_b;
  real_T c2_dm_y;
  real_T c2_ol_a;
  real_T c2_pl_b;
  real_T c2_em_y;
  real_T c2_av_x;
  real_T c2_bv_x;
  real_T c2_pl_a;
  real_T c2_ql_b;
  real_T c2_fm_y;
  real_T c2_cv_x;
  real_T c2_dv_x;
  real_T c2_ql_a;
  real_T c2_rl_b;
  real_T c2_gm_y;
  real_T c2_rl_a;
  real_T c2_sl_b;
  real_T c2_hm_y;
  real_T c2_ev_x;
  real_T c2_fv_x;
  real_T c2_sl_a;
  real_T c2_tl_b;
  real_T c2_im_y;
  real_T c2_gv_x;
  real_T c2_hv_x;
  real_T c2_tl_a;
  real_T c2_ul_b;
  real_T c2_jm_y;
  real_T c2_iv_x;
  real_T c2_jv_x;
  real_T c2_ul_a;
  real_T c2_vl_b;
  real_T c2_km_y;
  real_T c2_vl_a;
  real_T c2_wl_b;
  real_T c2_lm_y;
  real_T c2_kv_x;
  real_T c2_lv_x;
  real_T c2_wl_a;
  real_T c2_xl_b;
  real_T c2_mm_y;
  real_T c2_mv_x;
  real_T c2_nv_x;
  real_T c2_xl_a;
  real_T c2_yl_b;
  real_T c2_nm_y;
  real_T c2_yl_a;
  real_T c2_am_b;
  real_T c2_om_y;
  real_T c2_am_a;
  real_T c2_bm_b;
  real_T c2_pm_y;
  real_T c2_ov_x;
  real_T c2_pv_x;
  real_T c2_bm_a;
  real_T c2_cm_b;
  real_T c2_qm_y;
  real_T c2_qv_x;
  real_T c2_rv_x;
  real_T c2_cm_a;
  real_T c2_dm_b;
  real_T c2_rm_y;
  real_T c2_dm_a;
  real_T c2_em_b;
  real_T c2_sm_y;
  real_T c2_em_a;
  real_T c2_fm_b;
  real_T c2_tm_y;
  real_T c2_sv_x;
  real_T c2_tv_x;
  real_T c2_fm_a;
  real_T c2_gm_b;
  real_T c2_um_y;
  real_T c2_uv_x;
  real_T c2_vv_x;
  real_T c2_gm_a;
  real_T c2_hm_b;
  real_T c2_vm_y;
  real_T c2_hm_a;
  real_T c2_im_b;
  real_T c2_wm_y;
  real_T c2_wv_x;
  real_T c2_xv_x;
  real_T c2_im_a;
  real_T c2_jm_b;
  real_T c2_xm_y;
  real_T c2_jm_a;
  real_T c2_km_b;
  real_T c2_ym_y;
  real_T c2_km_a;
  real_T c2_lm_b;
  real_T c2_an_y;
  int32_T c2_i54;
  real_T c2_b_multipliers[9];
  real_T c2_b_delta1[3];
  real_T c2_dv0[3];
  int32_T c2_i55;
  real_T c2_g_A;
  real_T c2_B;
  real_T c2_yv_x;
  real_T c2_bn_y;
  real_T c2_aw_x;
  real_T c2_cn_y;
  real_T c2_h_A;
  real_T c2_b_B;
  real_T c2_bw_x;
  real_T c2_dn_y;
  real_T c2_cw_x;
  real_T c2_en_y;
  real_T c2_i_A;
  real_T c2_c_B;
  real_T c2_dw_x;
  real_T c2_fn_y;
  real_T c2_ew_x;
  real_T c2_gn_y;
  real_T c2_hn_y;
  real_T *c2_b_x_ddot;
  real_T *c2_b_y_ddot;
  real_T *c2_b_z_ddot;
  real_T *c2_b_phi_ddot;
  real_T *c2_b_theta_ddot;
  real_T *c2_b_psi_ddot;
  real_T (*c2_b_w)[4];
  real_T (*c2_b_params)[14];
  real_T (*c2_b_states)[12];
  c2_b_w = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
  c2_b_params = (real_T (*)[14])ssGetInputPortSignal(chartInstance->S, 1);
  c2_b_psi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c2_b_theta_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c2_b_phi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_b_z_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_b_y_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_b_x_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_b_states = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c2_sfEvent);
  for (c2_i3 = 0; c2_i3 < 12; c2_i3++) {
    c2_states[c2_i3] = (*c2_b_states)[c2_i3];
  }

  for (c2_i4 = 0; c2_i4 < 14; c2_i4++) {
    c2_params[c2_i4] = (*c2_b_params)[c2_i4];
  }

  for (c2_i5 = 0; c2_i5 < 4; c2_i5++) {
    c2_w[c2_i5] = (*c2_b_w)[c2_i5];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 61U, 61U, c2_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_m, 0U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_g, 1U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_l, 2U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Jxx, 3U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Jyy, 4U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Jzz, 5U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_k1, 6U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_k2, 7U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_k3, 8U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_k4, 9U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b1, 10U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b2, 11U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b3, 12U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b4, 13U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_w1, 14U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_w2, 15U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_w3, 16U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_w4, 17U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_x, 18U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_x_dot, 19U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_y, 20U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_y_dot, 21U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_z, 22U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_z_dot, 23U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_phi, 24U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_phi_dot, 25U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_theta, 26U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_theta_dot, 27U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_psi, 28U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_psi_dot, 29U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_multipliers, 30U, c2_f_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_upphi_x, 31U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_upphi_y, 32U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_upphi_z, 33U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_T_e2E, 34U, c2_f_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_force_mult_1, 35U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_force_mult_2, 36U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dx1_dphi, 37U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dx2_dphi, 38U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_upphi_phi, 39U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dx1_dtheta, 40U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dx2_dtheta, 41U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_upphi_theta, 42U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dx1_dpsi, 43U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dx2_dpsi, 44U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_upphi_psi, 45U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_delta1, 46U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_delta2, 47U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_delta3, 48U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_ddot, 49U, c2_e_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 50U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 51U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_states, 52U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_params, 53U, c2_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_w, 54U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_x_ddot, 55U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_y_ddot, 56U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_z_ddot, 57U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_phi_ddot, 58U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_theta_ddot, 59U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_psi_ddot, 60U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 3);
  c2_m = c2_params[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 4);
  c2_g = c2_params[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 5);
  c2_l = c2_params[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 7);
  c2_Jxx = c2_params[3];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 8);
  c2_Jyy = c2_params[4];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 9);
  c2_Jzz = c2_params[5];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 11);
  c2_k1 = c2_params[6];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 12);
  c2_k2 = c2_params[7];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 13);
  c2_k3 = c2_params[8];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 14);
  c2_k4 = c2_params[9];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 16);
  c2_b1 = c2_params[10];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 17);
  c2_b2 = c2_params[11];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 18);
  c2_b3 = c2_params[12];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 19);
  c2_b4 = c2_params[13];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 21);
  c2_w1 = c2_w[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 22);
  c2_w2 = c2_w[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 23);
  c2_w3 = c2_w[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 24);
  c2_w4 = c2_w[3];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 26);
  c2_x = c2_states[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 26);
  c2_x_dot = c2_states[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 26);
  c2_y = c2_states[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 26);
  c2_y_dot = c2_states[3];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 26);
  c2_z = c2_states[4];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 26);
  c2_z_dot = c2_states[5];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  c2_phi = c2_states[6];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  c2_phi_dot = c2_states[7];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  c2_theta = c2_states[8];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  c2_theta_dot = c2_states[9];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  c2_psi = c2_states[10];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  c2_psi_dot = c2_states[11];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 29);
  c2_b_x = c2_psi;
  c2_c_x = c2_b_x;
  c2_c_x = muDoubleScalarSin(c2_c_x);
  c2_a = c2_Jxx;
  c2_b = c2_mpower(chartInstance, c2_c_x);
  c2_b_y = c2_a * c2_b;
  c2_d_x = c2_theta;
  c2_e_x = c2_d_x;
  c2_e_x = muDoubleScalarSin(c2_e_x);
  c2_b_a = c2_b_y;
  c2_b_b = c2_mpower(chartInstance, c2_e_x);
  c2_c_y = c2_b_a * c2_b_b;
  c2_f_x = c2_psi;
  c2_g_x = c2_f_x;
  c2_g_x = muDoubleScalarCos(c2_g_x);
  c2_c_a = c2_Jyy;
  c2_c_b = c2_mpower(chartInstance, c2_g_x);
  c2_d_y = c2_c_a * c2_c_b;
  c2_h_x = c2_theta;
  c2_i_x = c2_h_x;
  c2_i_x = muDoubleScalarSin(c2_i_x);
  c2_d_a = c2_d_y;
  c2_d_b = c2_mpower(chartInstance, c2_i_x);
  c2_e_y = c2_d_a * c2_d_b;
  c2_j_x = c2_theta;
  c2_k_x = c2_j_x;
  c2_k_x = muDoubleScalarCos(c2_k_x);
  c2_e_a = c2_Jzz;
  c2_e_b = c2_mpower(chartInstance, c2_k_x);
  c2_f_y = c2_e_a * c2_e_b;
  c2_l_x = c2_psi;
  c2_m_x = c2_l_x;
  c2_m_x = muDoubleScalarSin(c2_m_x);
  c2_f_a = c2_Jxx;
  c2_f_b = c2_m_x;
  c2_g_y = c2_f_a * c2_f_b;
  c2_n_x = c2_theta;
  c2_o_x = c2_n_x;
  c2_o_x = muDoubleScalarSin(c2_o_x);
  c2_g_a = c2_g_y;
  c2_g_b = c2_o_x;
  c2_h_y = c2_g_a * c2_g_b;
  c2_p_x = c2_psi;
  c2_q_x = c2_p_x;
  c2_q_x = muDoubleScalarCos(c2_q_x);
  c2_h_a = c2_h_y;
  c2_h_b = c2_q_x;
  c2_i_y = c2_h_a * c2_h_b;
  c2_r_x = c2_psi;
  c2_s_x = c2_r_x;
  c2_s_x = muDoubleScalarCos(c2_s_x);
  c2_i_a = c2_Jyy;
  c2_i_b = c2_s_x;
  c2_j_y = c2_i_a * c2_i_b;
  c2_t_x = c2_theta;
  c2_u_x = c2_t_x;
  c2_u_x = muDoubleScalarSin(c2_u_x);
  c2_j_a = c2_j_y;
  c2_j_b = c2_u_x;
  c2_k_y = c2_j_a * c2_j_b;
  c2_v_x = c2_psi;
  c2_w_x = c2_v_x;
  c2_w_x = muDoubleScalarSin(c2_w_x);
  c2_k_a = c2_k_y;
  c2_k_b = c2_w_x;
  c2_l_y = c2_k_a * c2_k_b;
  c2_x_x = c2_theta;
  c2_y_x = c2_x_x;
  c2_y_x = muDoubleScalarCos(c2_y_x);
  c2_l_a = c2_Jzz;
  c2_l_b = c2_y_x;
  c2_m_y = c2_l_a * c2_l_b;
  c2_ab_x = c2_psi;
  c2_bb_x = c2_ab_x;
  c2_bb_x = muDoubleScalarSin(c2_bb_x);
  c2_m_a = c2_Jxx;
  c2_m_b = c2_bb_x;
  c2_n_y = c2_m_a * c2_m_b;
  c2_cb_x = c2_theta;
  c2_db_x = c2_cb_x;
  c2_db_x = muDoubleScalarSin(c2_db_x);
  c2_n_a = c2_n_y;
  c2_n_b = c2_db_x;
  c2_o_y = c2_n_a * c2_n_b;
  c2_eb_x = c2_psi;
  c2_fb_x = c2_eb_x;
  c2_fb_x = muDoubleScalarCos(c2_fb_x);
  c2_o_a = c2_o_y;
  c2_o_b = c2_fb_x;
  c2_p_y = c2_o_a * c2_o_b;
  c2_gb_x = c2_psi;
  c2_hb_x = c2_gb_x;
  c2_hb_x = muDoubleScalarCos(c2_hb_x);
  c2_p_a = c2_Jyy;
  c2_p_b = c2_hb_x;
  c2_q_y = c2_p_a * c2_p_b;
  c2_ib_x = c2_theta;
  c2_jb_x = c2_ib_x;
  c2_jb_x = muDoubleScalarSin(c2_jb_x);
  c2_q_a = c2_q_y;
  c2_q_b = c2_jb_x;
  c2_r_y = c2_q_a * c2_q_b;
  c2_kb_x = c2_psi;
  c2_lb_x = c2_kb_x;
  c2_lb_x = muDoubleScalarSin(c2_lb_x);
  c2_r_a = c2_r_y;
  c2_r_b = c2_lb_x;
  c2_s_y = c2_r_a * c2_r_b;
  c2_mb_x = c2_psi;
  c2_nb_x = c2_mb_x;
  c2_nb_x = muDoubleScalarCos(c2_nb_x);
  c2_s_a = c2_Jxx;
  c2_s_b = c2_mpower(chartInstance, c2_nb_x);
  c2_t_y = c2_s_a * c2_s_b;
  c2_ob_x = c2_psi;
  c2_pb_x = c2_ob_x;
  c2_pb_x = muDoubleScalarSin(c2_pb_x);
  c2_t_a = c2_Jyy;
  c2_t_b = c2_mpower(chartInstance, c2_pb_x);
  c2_u_y = c2_t_a * c2_t_b;
  c2_qb_x = c2_theta;
  c2_rb_x = c2_qb_x;
  c2_rb_x = muDoubleScalarCos(c2_rb_x);
  c2_u_a = c2_Jzz;
  c2_u_b = c2_rb_x;
  c2_v_y = c2_u_a * c2_u_b;
  c2_multipliers[0] = (c2_c_y + c2_e_y) + c2_f_y;
  c2_multipliers[3] = c2_i_y - c2_l_y;
  c2_multipliers[6] = c2_m_y;
  c2_multipliers[1] = c2_p_y - c2_s_y;
  c2_multipliers[4] = c2_t_y + c2_u_y;
  c2_multipliers[7] = 0.0;
  c2_multipliers[2] = c2_v_y;
  c2_multipliers[5] = 0.0;
  c2_multipliers[8] = c2_Jzz;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 35);
  c2_v_a = c2_k1;
  c2_v_b = c2_mpower(chartInstance, c2_w1);
  c2_w_y = c2_v_a * c2_v_b;
  c2_w_a = c2_k2;
  c2_w_b = c2_mpower(chartInstance, c2_w2);
  c2_x_y = c2_w_a * c2_w_b;
  c2_x_a = c2_k3;
  c2_x_b = c2_mpower(chartInstance, c2_w3);
  c2_y_y = c2_x_a * c2_x_b;
  c2_y_a = c2_k4;
  c2_y_b = c2_mpower(chartInstance, c2_w4);
  c2_ab_y = c2_y_a * c2_y_b;
  c2_sb_x = c2_theta;
  c2_tb_x = c2_sb_x;
  c2_tb_x = muDoubleScalarSin(c2_tb_x);
  c2_ab_a = ((c2_w_y + c2_x_y) + c2_y_y) + c2_ab_y;
  c2_ab_b = c2_tb_x;
  c2_bb_y = c2_ab_a * c2_ab_b;
  c2_ub_x = c2_phi;
  c2_vb_x = c2_ub_x;
  c2_vb_x = muDoubleScalarSin(c2_vb_x);
  c2_bb_a = c2_bb_y;
  c2_bb_b = c2_vb_x;
  c2_cb_y = c2_bb_a * c2_bb_b;
  c2_cb_a = -c2_b1;
  c2_cb_b = c2_mpower(chartInstance, c2_w1);
  c2_db_y = c2_cb_a * c2_cb_b;
  c2_db_a = c2_b2;
  c2_db_b = c2_mpower(chartInstance, c2_w2);
  c2_eb_y = c2_db_a * c2_db_b;
  c2_eb_a = c2_b3;
  c2_eb_b = c2_mpower(chartInstance, c2_w3);
  c2_fb_y = c2_eb_a * c2_eb_b;
  c2_fb_a = c2_b4;
  c2_fb_b = c2_mpower(chartInstance, c2_w4);
  c2_gb_y = c2_fb_a * c2_fb_b;
  c2_wb_x = c2_psi;
  c2_xb_x = c2_wb_x;
  c2_xb_x = muDoubleScalarCos(c2_xb_x);
  c2_yb_x = c2_phi;
  c2_ac_x = c2_yb_x;
  c2_ac_x = muDoubleScalarCos(c2_ac_x);
  c2_gb_a = c2_xb_x;
  c2_gb_b = c2_ac_x;
  c2_hb_y = c2_gb_a * c2_gb_b;
  c2_bc_x = c2_psi;
  c2_cc_x = c2_bc_x;
  c2_cc_x = muDoubleScalarSin(c2_cc_x);
  c2_dc_x = c2_theta;
  c2_ec_x = c2_dc_x;
  c2_ec_x = muDoubleScalarCos(c2_ec_x);
  c2_hb_a = c2_cc_x;
  c2_hb_b = c2_ec_x;
  c2_ib_y = c2_hb_a * c2_hb_b;
  c2_fc_x = c2_phi;
  c2_gc_x = c2_fc_x;
  c2_gc_x = muDoubleScalarSin(c2_gc_x);
  c2_ib_a = c2_ib_y;
  c2_ib_b = c2_gc_x;
  c2_jb_y = c2_ib_a * c2_ib_b;
  c2_jb_a = ((c2_db_y + c2_eb_y) + c2_fb_y) - c2_gb_y;
  c2_jb_b = c2_hb_y - c2_jb_y;
  c2_kb_y = c2_jb_a * c2_jb_b;
  c2_kb_a = c2_b1;
  c2_kb_b = c2_mpower(chartInstance, c2_w1);
  c2_lb_y = c2_kb_a * c2_kb_b;
  c2_lb_a = c2_b2;
  c2_lb_b = c2_mpower(chartInstance, c2_w2);
  c2_mb_y = c2_lb_a * c2_lb_b;
  c2_mb_a = c2_b3;
  c2_mb_b = c2_mpower(chartInstance, c2_w3);
  c2_nb_y = c2_mb_a * c2_mb_b;
  c2_nb_a = c2_b4;
  c2_nb_b = c2_mpower(chartInstance, c2_w4);
  c2_ob_y = c2_nb_a * c2_nb_b;
  c2_hc_x = c2_psi;
  c2_ic_x = c2_hc_x;
  c2_ic_x = muDoubleScalarSin(c2_ic_x);
  c2_jc_x = c2_phi;
  c2_kc_x = c2_jc_x;
  c2_kc_x = muDoubleScalarCos(c2_kc_x);
  c2_ob_a = -c2_ic_x;
  c2_ob_b = c2_kc_x;
  c2_pb_y = c2_ob_a * c2_ob_b;
  c2_lc_x = c2_psi;
  c2_mc_x = c2_lc_x;
  c2_mc_x = muDoubleScalarCos(c2_mc_x);
  c2_nc_x = c2_theta;
  c2_oc_x = c2_nc_x;
  c2_oc_x = muDoubleScalarCos(c2_oc_x);
  c2_pb_a = c2_mc_x;
  c2_pb_b = c2_oc_x;
  c2_qb_y = c2_pb_a * c2_pb_b;
  c2_pc_x = c2_phi;
  c2_qc_x = c2_pc_x;
  c2_qc_x = muDoubleScalarSin(c2_qc_x);
  c2_qb_a = c2_qb_y;
  c2_qb_b = c2_qc_x;
  c2_rb_y = c2_qb_a * c2_qb_b;
  c2_rb_a = ((c2_lb_y + c2_mb_y) - c2_nb_y) - c2_ob_y;
  c2_rb_b = c2_pb_y - c2_rb_y;
  c2_sb_y = c2_rb_a * c2_rb_b;
  c2_upphi_x = (c2_cb_y + c2_kb_y) + c2_sb_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 39);
  c2_sb_a = c2_k1;
  c2_sb_b = c2_mpower(chartInstance, c2_w1);
  c2_tb_y = c2_sb_a * c2_sb_b;
  c2_tb_a = c2_k2;
  c2_tb_b = c2_mpower(chartInstance, c2_w2);
  c2_ub_y = c2_tb_a * c2_tb_b;
  c2_ub_a = c2_k3;
  c2_ub_b = c2_mpower(chartInstance, c2_w3);
  c2_vb_y = c2_ub_a * c2_ub_b;
  c2_vb_a = c2_k4;
  c2_vb_b = c2_mpower(chartInstance, c2_w4);
  c2_wb_y = c2_vb_a * c2_vb_b;
  c2_rc_x = c2_theta;
  c2_sc_x = c2_rc_x;
  c2_sc_x = muDoubleScalarSin(c2_sc_x);
  c2_wb_a = -(((c2_tb_y + c2_ub_y) + c2_vb_y) + c2_wb_y);
  c2_wb_b = c2_sc_x;
  c2_xb_y = c2_wb_a * c2_wb_b;
  c2_tc_x = c2_phi;
  c2_uc_x = c2_tc_x;
  c2_uc_x = muDoubleScalarCos(c2_uc_x);
  c2_xb_a = c2_xb_y;
  c2_xb_b = c2_uc_x;
  c2_yb_y = c2_xb_a * c2_xb_b;
  c2_yb_a = -c2_b1;
  c2_yb_b = c2_mpower(chartInstance, c2_w1);
  c2_ac_y = c2_yb_a * c2_yb_b;
  c2_ac_a = c2_b2;
  c2_ac_b = c2_mpower(chartInstance, c2_w2);
  c2_bc_y = c2_ac_a * c2_ac_b;
  c2_bc_a = c2_b3;
  c2_bc_b = c2_mpower(chartInstance, c2_w3);
  c2_cc_y = c2_bc_a * c2_bc_b;
  c2_cc_a = c2_b4;
  c2_cc_b = c2_mpower(chartInstance, c2_w4);
  c2_dc_y = c2_cc_a * c2_cc_b;
  c2_vc_x = c2_psi;
  c2_wc_x = c2_vc_x;
  c2_wc_x = muDoubleScalarCos(c2_wc_x);
  c2_xc_x = c2_phi;
  c2_yc_x = c2_xc_x;
  c2_yc_x = muDoubleScalarSin(c2_yc_x);
  c2_dc_a = c2_wc_x;
  c2_dc_b = c2_yc_x;
  c2_ec_y = c2_dc_a * c2_dc_b;
  c2_ad_x = c2_psi;
  c2_bd_x = c2_ad_x;
  c2_bd_x = muDoubleScalarSin(c2_bd_x);
  c2_cd_x = c2_theta;
  c2_dd_x = c2_cd_x;
  c2_dd_x = muDoubleScalarCos(c2_dd_x);
  c2_ec_a = c2_bd_x;
  c2_ec_b = c2_dd_x;
  c2_fc_y = c2_ec_a * c2_ec_b;
  c2_ed_x = c2_phi;
  c2_fd_x = c2_ed_x;
  c2_fd_x = muDoubleScalarCos(c2_fd_x);
  c2_fc_a = c2_fc_y;
  c2_fc_b = c2_fd_x;
  c2_gc_y = c2_fc_a * c2_fc_b;
  c2_gc_a = ((c2_ac_y + c2_bc_y) + c2_cc_y) - c2_dc_y;
  c2_gc_b = c2_ec_y + c2_gc_y;
  c2_hc_y = c2_gc_a * c2_gc_b;
  c2_hc_a = c2_b1;
  c2_hc_b = c2_mpower(chartInstance, c2_w1);
  c2_ic_y = c2_hc_a * c2_hc_b;
  c2_ic_a = c2_b2;
  c2_ic_b = c2_mpower(chartInstance, c2_w2);
  c2_jc_y = c2_ic_a * c2_ic_b;
  c2_jc_a = c2_b3;
  c2_jc_b = c2_mpower(chartInstance, c2_w3);
  c2_kc_y = c2_jc_a * c2_jc_b;
  c2_kc_a = c2_b4;
  c2_kc_b = c2_mpower(chartInstance, c2_w4);
  c2_lc_y = c2_kc_a * c2_kc_b;
  c2_gd_x = c2_psi;
  c2_hd_x = c2_gd_x;
  c2_hd_x = muDoubleScalarSin(c2_hd_x);
  c2_id_x = c2_phi;
  c2_jd_x = c2_id_x;
  c2_jd_x = muDoubleScalarSin(c2_jd_x);
  c2_lc_a = -c2_hd_x;
  c2_lc_b = c2_jd_x;
  c2_mc_y = c2_lc_a * c2_lc_b;
  c2_kd_x = c2_psi;
  c2_ld_x = c2_kd_x;
  c2_ld_x = muDoubleScalarCos(c2_ld_x);
  c2_md_x = c2_theta;
  c2_nd_x = c2_md_x;
  c2_nd_x = muDoubleScalarCos(c2_nd_x);
  c2_mc_a = c2_ld_x;
  c2_mc_b = c2_nd_x;
  c2_nc_y = c2_mc_a * c2_mc_b;
  c2_od_x = c2_phi;
  c2_pd_x = c2_od_x;
  c2_pd_x = muDoubleScalarCos(c2_pd_x);
  c2_nc_a = c2_nc_y;
  c2_nc_b = c2_pd_x;
  c2_oc_y = c2_nc_a * c2_nc_b;
  c2_oc_a = ((c2_ic_y + c2_jc_y) - c2_kc_y) - c2_lc_y;
  c2_oc_b = c2_mc_y + c2_oc_y;
  c2_pc_y = c2_oc_a * c2_oc_b;
  c2_upphi_y = (c2_yb_y + c2_hc_y) + c2_pc_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 43);
  c2_pc_a = c2_k1;
  c2_pc_b = c2_mpower(chartInstance, c2_w1);
  c2_qc_y = c2_pc_a * c2_pc_b;
  c2_qc_a = c2_k2;
  c2_qc_b = c2_mpower(chartInstance, c2_w2);
  c2_rc_y = c2_qc_a * c2_qc_b;
  c2_rc_a = c2_k3;
  c2_rc_b = c2_mpower(chartInstance, c2_w3);
  c2_sc_y = c2_rc_a * c2_rc_b;
  c2_sc_a = c2_k4;
  c2_sc_b = c2_mpower(chartInstance, c2_w4);
  c2_tc_y = c2_sc_a * c2_sc_b;
  c2_qd_x = c2_theta;
  c2_rd_x = c2_qd_x;
  c2_rd_x = muDoubleScalarCos(c2_rd_x);
  c2_tc_a = ((c2_qc_y + c2_rc_y) + c2_sc_y) + c2_tc_y;
  c2_tc_b = c2_rd_x;
  c2_uc_y = c2_tc_a * c2_tc_b;
  c2_uc_a = -c2_b1;
  c2_uc_b = c2_mpower(chartInstance, c2_w1);
  c2_vc_y = c2_uc_a * c2_uc_b;
  c2_vc_a = c2_b2;
  c2_vc_b = c2_mpower(chartInstance, c2_w2);
  c2_wc_y = c2_vc_a * c2_vc_b;
  c2_wc_a = c2_b3;
  c2_wc_b = c2_mpower(chartInstance, c2_w3);
  c2_xc_y = c2_wc_a * c2_wc_b;
  c2_xc_a = c2_b4;
  c2_xc_b = c2_mpower(chartInstance, c2_w4);
  c2_yc_y = c2_xc_a * c2_xc_b;
  c2_sd_x = c2_psi;
  c2_td_x = c2_sd_x;
  c2_td_x = muDoubleScalarSin(c2_td_x);
  c2_yc_a = ((c2_vc_y + c2_wc_y) + c2_xc_y) - c2_yc_y;
  c2_yc_b = c2_td_x;
  c2_ad_y = c2_yc_a * c2_yc_b;
  c2_ud_x = c2_theta;
  c2_vd_x = c2_ud_x;
  c2_vd_x = muDoubleScalarSin(c2_vd_x);
  c2_ad_a = c2_ad_y;
  c2_ad_b = c2_vd_x;
  c2_bd_y = c2_ad_a * c2_ad_b;
  c2_bd_a = c2_b1;
  c2_bd_b = c2_mpower(chartInstance, c2_w1);
  c2_cd_y = c2_bd_a * c2_bd_b;
  c2_cd_a = c2_b2;
  c2_cd_b = c2_mpower(chartInstance, c2_w2);
  c2_dd_y = c2_cd_a * c2_cd_b;
  c2_dd_a = c2_b3;
  c2_dd_b = c2_mpower(chartInstance, c2_w3);
  c2_ed_y = c2_dd_a * c2_dd_b;
  c2_ed_a = c2_b4;
  c2_ed_b = c2_mpower(chartInstance, c2_w4);
  c2_fd_y = c2_ed_a * c2_ed_b;
  c2_wd_x = c2_psi;
  c2_xd_x = c2_wd_x;
  c2_xd_x = muDoubleScalarCos(c2_xd_x);
  c2_fd_a = ((c2_cd_y + c2_dd_y) - c2_ed_y) - c2_fd_y;
  c2_fd_b = c2_xd_x;
  c2_gd_y = c2_fd_a * c2_fd_b;
  c2_yd_x = c2_theta;
  c2_ae_x = c2_yd_x;
  c2_ae_x = muDoubleScalarSin(c2_ae_x);
  c2_gd_a = c2_gd_y;
  c2_gd_b = c2_ae_x;
  c2_hd_y = c2_gd_a * c2_gd_b;
  c2_upphi_z = (c2_uc_y + c2_bd_y) + c2_hd_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 47);
  c2_be_x = c2_psi;
  c2_ce_x = c2_be_x;
  c2_ce_x = muDoubleScalarCos(c2_ce_x);
  c2_de_x = c2_phi;
  c2_ee_x = c2_de_x;
  c2_ee_x = muDoubleScalarCos(c2_ee_x);
  c2_hd_a = c2_ce_x;
  c2_hd_b = c2_ee_x;
  c2_id_y = c2_hd_a * c2_hd_b;
  c2_fe_x = c2_psi;
  c2_ge_x = c2_fe_x;
  c2_ge_x = muDoubleScalarSin(c2_ge_x);
  c2_he_x = c2_theta;
  c2_ie_x = c2_he_x;
  c2_ie_x = muDoubleScalarCos(c2_ie_x);
  c2_id_a = c2_ge_x;
  c2_id_b = c2_ie_x;
  c2_jd_y = c2_id_a * c2_id_b;
  c2_je_x = c2_phi;
  c2_ke_x = c2_je_x;
  c2_ke_x = muDoubleScalarSin(c2_ke_x);
  c2_jd_a = c2_jd_y;
  c2_jd_b = c2_ke_x;
  c2_kd_y = c2_jd_a * c2_jd_b;
  c2_le_x = c2_psi;
  c2_me_x = c2_le_x;
  c2_me_x = muDoubleScalarSin(c2_me_x);
  c2_ne_x = c2_phi;
  c2_oe_x = c2_ne_x;
  c2_oe_x = muDoubleScalarCos(c2_oe_x);
  c2_kd_a = -c2_me_x;
  c2_kd_b = c2_oe_x;
  c2_ld_y = c2_kd_a * c2_kd_b;
  c2_pe_x = c2_psi;
  c2_qe_x = c2_pe_x;
  c2_qe_x = muDoubleScalarCos(c2_qe_x);
  c2_re_x = c2_theta;
  c2_se_x = c2_re_x;
  c2_se_x = muDoubleScalarCos(c2_se_x);
  c2_ld_a = c2_qe_x;
  c2_ld_b = c2_se_x;
  c2_md_y = c2_ld_a * c2_ld_b;
  c2_te_x = c2_phi;
  c2_ue_x = c2_te_x;
  c2_ue_x = muDoubleScalarSin(c2_ue_x);
  c2_md_a = c2_md_y;
  c2_md_b = c2_ue_x;
  c2_nd_y = c2_md_a * c2_md_b;
  c2_ve_x = c2_theta;
  c2_we_x = c2_ve_x;
  c2_we_x = muDoubleScalarSin(c2_we_x);
  c2_xe_x = c2_phi;
  c2_ye_x = c2_xe_x;
  c2_ye_x = muDoubleScalarSin(c2_ye_x);
  c2_nd_a = c2_we_x;
  c2_nd_b = c2_ye_x;
  c2_od_y = c2_nd_a * c2_nd_b;
  c2_af_x = c2_psi;
  c2_bf_x = c2_af_x;
  c2_bf_x = muDoubleScalarCos(c2_bf_x);
  c2_cf_x = c2_phi;
  c2_df_x = c2_cf_x;
  c2_df_x = muDoubleScalarSin(c2_df_x);
  c2_od_a = c2_bf_x;
  c2_od_b = c2_df_x;
  c2_pd_y = c2_od_a * c2_od_b;
  c2_ef_x = c2_psi;
  c2_ff_x = c2_ef_x;
  c2_ff_x = muDoubleScalarSin(c2_ff_x);
  c2_gf_x = c2_theta;
  c2_hf_x = c2_gf_x;
  c2_hf_x = muDoubleScalarCos(c2_hf_x);
  c2_pd_a = c2_ff_x;
  c2_pd_b = c2_hf_x;
  c2_qd_y = c2_pd_a * c2_pd_b;
  c2_if_x = c2_phi;
  c2_jf_x = c2_if_x;
  c2_jf_x = muDoubleScalarCos(c2_jf_x);
  c2_qd_a = c2_qd_y;
  c2_qd_b = c2_jf_x;
  c2_rd_y = c2_qd_a * c2_qd_b;
  c2_kf_x = c2_psi;
  c2_lf_x = c2_kf_x;
  c2_lf_x = muDoubleScalarSin(c2_lf_x);
  c2_mf_x = c2_phi;
  c2_nf_x = c2_mf_x;
  c2_nf_x = muDoubleScalarSin(c2_nf_x);
  c2_rd_a = -c2_lf_x;
  c2_rd_b = c2_nf_x;
  c2_sd_y = c2_rd_a * c2_rd_b;
  c2_of_x = c2_psi;
  c2_pf_x = c2_of_x;
  c2_pf_x = muDoubleScalarCos(c2_pf_x);
  c2_qf_x = c2_theta;
  c2_rf_x = c2_qf_x;
  c2_rf_x = muDoubleScalarCos(c2_rf_x);
  c2_sd_a = c2_pf_x;
  c2_sd_b = c2_rf_x;
  c2_td_y = c2_sd_a * c2_sd_b;
  c2_sf_x = c2_phi;
  c2_tf_x = c2_sf_x;
  c2_tf_x = muDoubleScalarCos(c2_tf_x);
  c2_td_a = c2_td_y;
  c2_td_b = c2_tf_x;
  c2_ud_y = c2_td_a * c2_td_b;
  c2_uf_x = c2_theta;
  c2_vf_x = c2_uf_x;
  c2_vf_x = muDoubleScalarSin(c2_vf_x);
  c2_wf_x = c2_phi;
  c2_xf_x = c2_wf_x;
  c2_xf_x = muDoubleScalarCos(c2_xf_x);
  c2_ud_a = -c2_vf_x;
  c2_ud_b = c2_xf_x;
  c2_vd_y = c2_ud_a * c2_ud_b;
  c2_yf_x = c2_psi;
  c2_ag_x = c2_yf_x;
  c2_ag_x = muDoubleScalarSin(c2_ag_x);
  c2_bg_x = c2_theta;
  c2_cg_x = c2_bg_x;
  c2_cg_x = muDoubleScalarSin(c2_cg_x);
  c2_vd_a = c2_ag_x;
  c2_vd_b = c2_cg_x;
  c2_wd_y = c2_vd_a * c2_vd_b;
  c2_dg_x = c2_psi;
  c2_eg_x = c2_dg_x;
  c2_eg_x = muDoubleScalarCos(c2_eg_x);
  c2_fg_x = c2_theta;
  c2_gg_x = c2_fg_x;
  c2_gg_x = muDoubleScalarSin(c2_gg_x);
  c2_wd_a = c2_eg_x;
  c2_wd_b = c2_gg_x;
  c2_xd_y = c2_wd_a * c2_wd_b;
  c2_hg_x = c2_theta;
  c2_ig_x = c2_hg_x;
  c2_ig_x = muDoubleScalarCos(c2_ig_x);
  c2_T_e2E[0] = c2_id_y - c2_kd_y;
  c2_T_e2E[3] = c2_ld_y - c2_nd_y;
  c2_T_e2E[6] = c2_od_y;
  c2_T_e2E[1] = c2_pd_y + c2_rd_y;
  c2_T_e2E[4] = c2_sd_y + c2_ud_y;
  c2_T_e2E[7] = c2_vd_y;
  c2_T_e2E[2] = c2_wd_y;
  c2_T_e2E[5] = c2_xd_y;
  c2_T_e2E[8] = c2_ig_x;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 52);
  c2_xd_a = -c2_b1;
  c2_xd_b = c2_mpower(chartInstance, c2_w1);
  c2_yd_y = c2_xd_a * c2_xd_b;
  c2_yd_a = c2_b3;
  c2_yd_b = c2_mpower(chartInstance, c2_w3);
  c2_ae_y = c2_yd_a * c2_yd_b;
  c2_ae_a = c2_b1;
  c2_ae_b = c2_mpower(chartInstance, c2_w1);
  c2_be_y = c2_ae_a * c2_ae_b;
  c2_be_a = c2_b3;
  c2_be_b = c2_mpower(chartInstance, c2_w3);
  c2_ce_y = c2_be_a * c2_be_b;
  c2_ce_a = c2_k1;
  c2_ce_b = c2_mpower(chartInstance, c2_w1);
  c2_de_y = c2_ce_a * c2_ce_b;
  c2_de_a = c2_k3;
  c2_de_b = c2_mpower(chartInstance, c2_w3);
  c2_ee_y = c2_de_a * c2_de_b;
  c2_force_mult_1[0] = c2_yd_y - c2_ae_y;
  c2_force_mult_1[1] = c2_be_y + c2_ce_y;
  c2_force_mult_1[2] = c2_de_y - c2_ee_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 53);
  c2_ee_a = c2_b2;
  c2_ee_b = c2_mpower(chartInstance, c2_w2);
  c2_fe_y = c2_ee_a * c2_ee_b;
  c2_fe_a = c2_b4;
  c2_fe_b = c2_mpower(chartInstance, c2_w4);
  c2_ge_y = c2_fe_a * c2_fe_b;
  c2_ge_a = c2_b2;
  c2_ge_b = c2_mpower(chartInstance, c2_w2);
  c2_he_y = c2_ge_a * c2_ge_b;
  c2_he_a = c2_b4;
  c2_he_b = c2_mpower(chartInstance, c2_w4);
  c2_ie_y = c2_he_a * c2_he_b;
  c2_ie_a = c2_k2;
  c2_ie_b = c2_mpower(chartInstance, c2_w2);
  c2_je_y = c2_ie_a * c2_ie_b;
  c2_je_a = c2_k4;
  c2_je_b = c2_mpower(chartInstance, c2_w4);
  c2_ke_y = c2_je_a * c2_je_b;
  c2_force_mult_2[0] = c2_fe_y + c2_ge_y;
  c2_force_mult_2[1] = c2_he_y + c2_ie_y;
  c2_force_mult_2[2] = c2_je_y - c2_ke_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 55);
  c2_A = c2_l;
  c2_jg_x = c2_A;
  c2_kg_x = c2_jg_x;
  c2_le_y = c2_kg_x / 1.4142135623730951;
  c2_lg_x = c2_psi;
  c2_mg_x = c2_lg_x;
  c2_mg_x = muDoubleScalarSin(c2_mg_x);
  c2_ng_x = c2_phi;
  c2_og_x = c2_ng_x;
  c2_og_x = muDoubleScalarSin(c2_og_x);
  c2_ke_a = c2_mg_x;
  c2_ke_b = c2_og_x;
  c2_me_y = c2_ke_a * c2_ke_b;
  c2_pg_x = c2_psi;
  c2_qg_x = c2_pg_x;
  c2_qg_x = muDoubleScalarCos(c2_qg_x);
  c2_rg_x = c2_theta;
  c2_sg_x = c2_rg_x;
  c2_sg_x = muDoubleScalarCos(c2_sg_x);
  c2_le_a = c2_qg_x;
  c2_le_b = c2_sg_x;
  c2_ne_y = c2_le_a * c2_le_b;
  c2_tg_x = c2_phi;
  c2_ug_x = c2_tg_x;
  c2_ug_x = muDoubleScalarCos(c2_ug_x);
  c2_me_a = c2_ne_y;
  c2_me_b = c2_ug_x;
  c2_oe_y = c2_me_a * c2_me_b;
  c2_vg_x = c2_psi;
  c2_wg_x = c2_vg_x;
  c2_wg_x = muDoubleScalarCos(c2_wg_x);
  c2_xg_x = c2_phi;
  c2_yg_x = c2_xg_x;
  c2_yg_x = muDoubleScalarSin(c2_yg_x);
  c2_ne_a = c2_wg_x;
  c2_ne_b = c2_yg_x;
  c2_pe_y = c2_ne_a * c2_ne_b;
  c2_ah_x = c2_psi;
  c2_bh_x = c2_ah_x;
  c2_bh_x = muDoubleScalarSin(c2_bh_x);
  c2_ch_x = c2_theta;
  c2_dh_x = c2_ch_x;
  c2_dh_x = muDoubleScalarCos(c2_dh_x);
  c2_oe_a = c2_bh_x;
  c2_oe_b = c2_dh_x;
  c2_qe_y = c2_oe_a * c2_oe_b;
  c2_eh_x = c2_phi;
  c2_fh_x = c2_eh_x;
  c2_fh_x = muDoubleScalarCos(c2_fh_x);
  c2_pe_a = c2_qe_y;
  c2_pe_b = c2_fh_x;
  c2_re_y = c2_pe_a * c2_pe_b;
  c2_gh_x = c2_psi;
  c2_hh_x = c2_gh_x;
  c2_hh_x = muDoubleScalarCos(c2_hh_x);
  c2_ih_x = c2_phi;
  c2_jh_x = c2_ih_x;
  c2_jh_x = muDoubleScalarCos(c2_jh_x);
  c2_qe_a = c2_hh_x;
  c2_qe_b = c2_jh_x;
  c2_se_y = c2_qe_a * c2_qe_b;
  c2_kh_x = c2_psi;
  c2_lh_x = c2_kh_x;
  c2_lh_x = muDoubleScalarSin(c2_lh_x);
  c2_mh_x = c2_theta;
  c2_nh_x = c2_mh_x;
  c2_nh_x = muDoubleScalarCos(c2_nh_x);
  c2_re_a = c2_lh_x;
  c2_re_b = c2_nh_x;
  c2_te_y = c2_re_a * c2_re_b;
  c2_oh_x = c2_phi;
  c2_ph_x = c2_oh_x;
  c2_ph_x = muDoubleScalarSin(c2_ph_x);
  c2_se_a = c2_te_y;
  c2_se_b = c2_ph_x;
  c2_ue_y = c2_se_a * c2_se_b;
  c2_qh_x = c2_psi;
  c2_rh_x = c2_qh_x;
  c2_rh_x = muDoubleScalarSin(c2_rh_x);
  c2_sh_x = c2_phi;
  c2_th_x = c2_sh_x;
  c2_th_x = muDoubleScalarCos(c2_th_x);
  c2_te_a = c2_rh_x;
  c2_te_b = c2_th_x;
  c2_ve_y = c2_te_a * c2_te_b;
  c2_uh_x = c2_psi;
  c2_vh_x = c2_uh_x;
  c2_vh_x = muDoubleScalarCos(c2_vh_x);
  c2_wh_x = c2_theta;
  c2_xh_x = c2_wh_x;
  c2_xh_x = muDoubleScalarCos(c2_xh_x);
  c2_ue_a = c2_vh_x;
  c2_ue_b = c2_xh_x;
  c2_we_y = c2_ue_a * c2_ue_b;
  c2_yh_x = c2_phi;
  c2_ai_x = c2_yh_x;
  c2_ai_x = muDoubleScalarSin(c2_ai_x);
  c2_ve_a = c2_we_y;
  c2_ve_b = c2_ai_x;
  c2_xe_y = c2_ve_a * c2_ve_b;
  c2_we_a = c2_le_y;
  c2_we_b[0] = ((c2_me_y - c2_oe_y) - c2_pe_y) - c2_re_y;
  c2_we_b[1] = ((c2_se_y - c2_ue_y) - c2_ve_y) - c2_xe_y;
  c2_we_b[2] = 0.0;
  for (c2_i6 = 0; c2_i6 < 3; c2_i6++) {
    c2_dx1_dphi[c2_i6] = c2_we_a * c2_we_b[c2_i6];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 60);
  c2_b_A = c2_l;
  c2_bi_x = c2_b_A;
  c2_ci_x = c2_bi_x;
  c2_ye_y = c2_ci_x / 1.4142135623730951;
  c2_di_x = c2_psi;
  c2_ei_x = c2_di_x;
  c2_ei_x = muDoubleScalarSin(c2_ei_x);
  c2_fi_x = c2_phi;
  c2_gi_x = c2_fi_x;
  c2_gi_x = muDoubleScalarSin(c2_gi_x);
  c2_xe_a = c2_ei_x;
  c2_xe_b = c2_gi_x;
  c2_af_y = c2_xe_a * c2_xe_b;
  c2_hi_x = c2_psi;
  c2_ii_x = c2_hi_x;
  c2_ii_x = muDoubleScalarCos(c2_ii_x);
  c2_ji_x = c2_theta;
  c2_ki_x = c2_ji_x;
  c2_ki_x = muDoubleScalarCos(c2_ki_x);
  c2_ye_a = c2_ii_x;
  c2_ye_b = c2_ki_x;
  c2_bf_y = c2_ye_a * c2_ye_b;
  c2_li_x = c2_phi;
  c2_mi_x = c2_li_x;
  c2_mi_x = muDoubleScalarCos(c2_mi_x);
  c2_af_a = c2_bf_y;
  c2_af_b = c2_mi_x;
  c2_cf_y = c2_af_a * c2_af_b;
  c2_ni_x = c2_psi;
  c2_oi_x = c2_ni_x;
  c2_oi_x = muDoubleScalarCos(c2_oi_x);
  c2_pi_x = c2_phi;
  c2_qi_x = c2_pi_x;
  c2_qi_x = muDoubleScalarSin(c2_qi_x);
  c2_bf_a = c2_oi_x;
  c2_bf_b = c2_qi_x;
  c2_df_y = c2_bf_a * c2_bf_b;
  c2_ri_x = c2_psi;
  c2_si_x = c2_ri_x;
  c2_si_x = muDoubleScalarSin(c2_si_x);
  c2_ti_x = c2_theta;
  c2_ui_x = c2_ti_x;
  c2_ui_x = muDoubleScalarCos(c2_ui_x);
  c2_cf_a = c2_si_x;
  c2_cf_b = c2_ui_x;
  c2_ef_y = c2_cf_a * c2_cf_b;
  c2_vi_x = c2_phi;
  c2_wi_x = c2_vi_x;
  c2_wi_x = muDoubleScalarCos(c2_wi_x);
  c2_df_a = c2_ef_y;
  c2_df_b = c2_wi_x;
  c2_ff_y = c2_df_a * c2_df_b;
  c2_xi_x = c2_psi;
  c2_yi_x = c2_xi_x;
  c2_yi_x = muDoubleScalarSin(c2_yi_x);
  c2_aj_x = c2_phi;
  c2_bj_x = c2_aj_x;
  c2_bj_x = muDoubleScalarCos(c2_bj_x);
  c2_ef_a = -c2_yi_x;
  c2_ef_b = c2_bj_x;
  c2_gf_y = c2_ef_a * c2_ef_b;
  c2_cj_x = c2_psi;
  c2_dj_x = c2_cj_x;
  c2_dj_x = muDoubleScalarCos(c2_dj_x);
  c2_ej_x = c2_theta;
  c2_fj_x = c2_ej_x;
  c2_fj_x = muDoubleScalarCos(c2_fj_x);
  c2_ff_a = c2_dj_x;
  c2_ff_b = c2_fj_x;
  c2_hf_y = c2_ff_a * c2_ff_b;
  c2_gj_x = c2_phi;
  c2_hj_x = c2_gj_x;
  c2_hj_x = muDoubleScalarSin(c2_hj_x);
  c2_gf_a = c2_hf_y;
  c2_gf_b = c2_hj_x;
  c2_if_y = c2_gf_a * c2_gf_b;
  c2_ij_x = c2_psi;
  c2_jj_x = c2_ij_x;
  c2_jj_x = muDoubleScalarCos(c2_jj_x);
  c2_kj_x = c2_phi;
  c2_lj_x = c2_kj_x;
  c2_lj_x = muDoubleScalarCos(c2_lj_x);
  c2_hf_a = c2_jj_x;
  c2_hf_b = c2_lj_x;
  c2_jf_y = c2_hf_a * c2_hf_b;
  c2_mj_x = c2_psi;
  c2_nj_x = c2_mj_x;
  c2_nj_x = muDoubleScalarSin(c2_nj_x);
  c2_oj_x = c2_theta;
  c2_pj_x = c2_oj_x;
  c2_pj_x = muDoubleScalarCos(c2_pj_x);
  c2_if_a = c2_nj_x;
  c2_if_b = c2_pj_x;
  c2_kf_y = c2_if_a * c2_if_b;
  c2_qj_x = c2_phi;
  c2_rj_x = c2_qj_x;
  c2_rj_x = muDoubleScalarSin(c2_rj_x);
  c2_jf_a = c2_kf_y;
  c2_jf_b = c2_rj_x;
  c2_lf_y = c2_jf_a * c2_jf_b;
  c2_kf_a = c2_ye_y;
  c2_we_b[0] = ((c2_af_y - c2_cf_y) + c2_df_y) + c2_ff_y;
  c2_we_b[1] = ((c2_gf_y - c2_if_y) - c2_jf_y) + c2_lf_y;
  c2_we_b[2] = 0.0;
  for (c2_i7 = 0; c2_i7 < 3; c2_i7++) {
    c2_dx2_dphi[c2_i7] = c2_kf_a * c2_we_b[c2_i7];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 65);
  for (c2_i8 = 0; c2_i8 < 9; c2_i8++) {
    c2_lf_a[c2_i8] = c2_T_e2E[c2_i8];
  }

  for (c2_i9 = 0; c2_i9 < 3; c2_i9++) {
    c2_we_b[c2_i9] = c2_force_mult_1[c2_i9];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i10 = 0; c2_i10 < 3; c2_i10++) {
    c2_mf_y[c2_i10] = 0.0;
    c2_i11 = 0;
    for (c2_i12 = 0; c2_i12 < 3; c2_i12++) {
      c2_mf_y[c2_i10] += c2_lf_a[c2_i11 + c2_i10] * c2_we_b[c2_i12];
      c2_i11 += 3;
    }
  }

  for (c2_i13 = 0; c2_i13 < 9; c2_i13++) {
    c2_lf_a[c2_i13] = c2_T_e2E[c2_i13];
  }

  for (c2_i14 = 0; c2_i14 < 3; c2_i14++) {
    c2_we_b[c2_i14] = c2_force_mult_2[c2_i14];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i15 = 0; c2_i15 < 3; c2_i15++) {
    c2_nf_y[c2_i15] = 0.0;
    c2_i16 = 0;
    for (c2_i17 = 0; c2_i17 < 3; c2_i17++) {
      c2_nf_y[c2_i15] += c2_lf_a[c2_i16 + c2_i15] * c2_we_b[c2_i17];
      c2_i16 += 3;
    }
  }

  for (c2_i18 = 0; c2_i18 < 3; c2_i18++) {
    c2_of_y[c2_i18] = c2_mf_y[c2_i18];
  }

  for (c2_i19 = 0; c2_i19 < 3; c2_i19++) {
    c2_b_dx1_dphi[c2_i19] = c2_dx1_dphi[c2_i19];
  }

  for (c2_i20 = 0; c2_i20 < 3; c2_i20++) {
    c2_pf_y[c2_i20] = c2_nf_y[c2_i20];
  }

  for (c2_i21 = 0; c2_i21 < 3; c2_i21++) {
    c2_b_dx2_dphi[c2_i21] = c2_dx2_dphi[c2_i21];
  }

  c2_upphi_phi = c2_dot(chartInstance, c2_of_y, c2_b_dx1_dphi) + c2_dot
    (chartInstance, c2_pf_y, c2_b_dx2_dphi);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 68);
  c2_c_A = c2_l;
  c2_sj_x = c2_c_A;
  c2_tj_x = c2_sj_x;
  c2_qf_y = c2_tj_x / 1.4142135623730951;
  c2_uj_x = c2_psi;
  c2_vj_x = c2_uj_x;
  c2_vj_x = muDoubleScalarSin(c2_vj_x);
  c2_wj_x = c2_theta;
  c2_xj_x = c2_wj_x;
  c2_xj_x = muDoubleScalarSin(c2_xj_x);
  c2_mf_a = c2_vj_x;
  c2_kf_b = c2_xj_x;
  c2_rf_y = c2_mf_a * c2_kf_b;
  c2_yj_x = c2_phi;
  c2_ak_x = c2_yj_x;
  c2_ak_x = muDoubleScalarSin(c2_ak_x);
  c2_nf_a = c2_rf_y;
  c2_lf_b = c2_ak_x;
  c2_sf_y = c2_nf_a * c2_lf_b;
  c2_bk_x = c2_psi;
  c2_ck_x = c2_bk_x;
  c2_ck_x = muDoubleScalarCos(c2_ck_x);
  c2_dk_x = c2_theta;
  c2_ek_x = c2_dk_x;
  c2_ek_x = muDoubleScalarSin(c2_ek_x);
  c2_of_a = c2_ck_x;
  c2_mf_b = c2_ek_x;
  c2_tf_y = c2_of_a * c2_mf_b;
  c2_fk_x = c2_phi;
  c2_gk_x = c2_fk_x;
  c2_gk_x = muDoubleScalarSin(c2_gk_x);
  c2_pf_a = c2_tf_y;
  c2_nf_b = c2_gk_x;
  c2_uf_y = c2_pf_a * c2_nf_b;
  c2_hk_x = c2_psi;
  c2_ik_x = c2_hk_x;
  c2_ik_x = muDoubleScalarSin(c2_ik_x);
  c2_jk_x = c2_theta;
  c2_kk_x = c2_jk_x;
  c2_kk_x = muDoubleScalarSin(c2_kk_x);
  c2_qf_a = -c2_ik_x;
  c2_of_b = c2_kk_x;
  c2_vf_y = c2_qf_a * c2_of_b;
  c2_lk_x = c2_phi;
  c2_mk_x = c2_lk_x;
  c2_mk_x = muDoubleScalarCos(c2_mk_x);
  c2_rf_a = c2_vf_y;
  c2_pf_b = c2_mk_x;
  c2_wf_y = c2_rf_a * c2_pf_b;
  c2_nk_x = c2_psi;
  c2_ok_x = c2_nk_x;
  c2_ok_x = muDoubleScalarCos(c2_ok_x);
  c2_pk_x = c2_theta;
  c2_qk_x = c2_pk_x;
  c2_qk_x = muDoubleScalarSin(c2_qk_x);
  c2_sf_a = c2_ok_x;
  c2_qf_b = c2_qk_x;
  c2_xf_y = c2_sf_a * c2_qf_b;
  c2_rk_x = c2_phi;
  c2_sk_x = c2_rk_x;
  c2_sk_x = muDoubleScalarCos(c2_sk_x);
  c2_tf_a = c2_xf_y;
  c2_rf_b = c2_sk_x;
  c2_yf_y = c2_tf_a * c2_rf_b;
  c2_tk_x = c2_psi;
  c2_uk_x = c2_tk_x;
  c2_uk_x = muDoubleScalarSin(c2_uk_x);
  c2_vk_x = c2_theta;
  c2_wk_x = c2_vk_x;
  c2_wk_x = muDoubleScalarCos(c2_wk_x);
  c2_uf_a = c2_uk_x;
  c2_sf_b = c2_wk_x;
  c2_ag_y = c2_uf_a * c2_sf_b;
  c2_xk_x = c2_psi;
  c2_yk_x = c2_xk_x;
  c2_yk_x = muDoubleScalarCos(c2_yk_x);
  c2_al_x = c2_theta;
  c2_bl_x = c2_al_x;
  c2_bl_x = muDoubleScalarCos(c2_bl_x);
  c2_vf_a = c2_yk_x;
  c2_tf_b = c2_bl_x;
  c2_bg_y = c2_vf_a * c2_tf_b;
  c2_wf_a = c2_qf_y;
  c2_we_b[0] = c2_sf_y + c2_uf_y;
  c2_we_b[1] = c2_wf_y - c2_yf_y;
  c2_we_b[2] = c2_ag_y + c2_bg_y;
  for (c2_i22 = 0; c2_i22 < 3; c2_i22++) {
    c2_dx1_dtheta[c2_i22] = c2_wf_a * c2_we_b[c2_i22];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 73);
  c2_d_A = c2_l;
  c2_cl_x = c2_d_A;
  c2_dl_x = c2_cl_x;
  c2_cg_y = c2_dl_x / 1.4142135623730951;
  c2_el_x = c2_psi;
  c2_fl_x = c2_el_x;
  c2_fl_x = muDoubleScalarCos(c2_fl_x);
  c2_gl_x = c2_theta;
  c2_hl_x = c2_gl_x;
  c2_hl_x = muDoubleScalarSin(c2_hl_x);
  c2_xf_a = c2_fl_x;
  c2_uf_b = c2_hl_x;
  c2_dg_y = c2_xf_a * c2_uf_b;
  c2_il_x = c2_phi;
  c2_jl_x = c2_il_x;
  c2_jl_x = muDoubleScalarSin(c2_jl_x);
  c2_yf_a = c2_dg_y;
  c2_vf_b = c2_jl_x;
  c2_eg_y = c2_yf_a * c2_vf_b;
  c2_kl_x = c2_psi;
  c2_ll_x = c2_kl_x;
  c2_ll_x = muDoubleScalarSin(c2_ll_x);
  c2_ml_x = c2_theta;
  c2_nl_x = c2_ml_x;
  c2_nl_x = muDoubleScalarSin(c2_nl_x);
  c2_ag_a = c2_ll_x;
  c2_wf_b = c2_nl_x;
  c2_fg_y = c2_ag_a * c2_wf_b;
  c2_ol_x = c2_phi;
  c2_pl_x = c2_ol_x;
  c2_pl_x = muDoubleScalarSin(c2_pl_x);
  c2_bg_a = c2_fg_y;
  c2_xf_b = c2_pl_x;
  c2_gg_y = c2_bg_a * c2_xf_b;
  c2_ql_x = c2_psi;
  c2_rl_x = c2_ql_x;
  c2_rl_x = muDoubleScalarCos(c2_rl_x);
  c2_sl_x = c2_theta;
  c2_tl_x = c2_sl_x;
  c2_tl_x = muDoubleScalarSin(c2_tl_x);
  c2_cg_a = -c2_rl_x;
  c2_yf_b = c2_tl_x;
  c2_hg_y = c2_cg_a * c2_yf_b;
  c2_ul_x = c2_phi;
  c2_vl_x = c2_ul_x;
  c2_vl_x = muDoubleScalarCos(c2_vl_x);
  c2_dg_a = c2_hg_y;
  c2_ag_b = c2_vl_x;
  c2_ig_y = c2_dg_a * c2_ag_b;
  c2_wl_x = c2_psi;
  c2_xl_x = c2_wl_x;
  c2_xl_x = muDoubleScalarSin(c2_xl_x);
  c2_yl_x = c2_theta;
  c2_am_x = c2_yl_x;
  c2_am_x = muDoubleScalarSin(c2_am_x);
  c2_eg_a = c2_xl_x;
  c2_bg_b = c2_am_x;
  c2_jg_y = c2_eg_a * c2_bg_b;
  c2_bm_x = c2_phi;
  c2_cm_x = c2_bm_x;
  c2_cm_x = muDoubleScalarCos(c2_cm_x);
  c2_fg_a = c2_jg_y;
  c2_cg_b = c2_cm_x;
  c2_kg_y = c2_fg_a * c2_cg_b;
  c2_dm_x = c2_psi;
  c2_em_x = c2_dm_x;
  c2_em_x = muDoubleScalarCos(c2_em_x);
  c2_fm_x = c2_theta;
  c2_gm_x = c2_fm_x;
  c2_gm_x = muDoubleScalarCos(c2_gm_x);
  c2_gg_a = c2_em_x;
  c2_dg_b = c2_gm_x;
  c2_lg_y = c2_gg_a * c2_dg_b;
  c2_hm_x = c2_psi;
  c2_im_x = c2_hm_x;
  c2_im_x = muDoubleScalarSin(c2_im_x);
  c2_jm_x = c2_theta;
  c2_km_x = c2_jm_x;
  c2_km_x = muDoubleScalarCos(c2_km_x);
  c2_hg_a = c2_im_x;
  c2_eg_b = c2_km_x;
  c2_mg_y = c2_hg_a * c2_eg_b;
  c2_ig_a = c2_cg_y;
  c2_we_b[0] = c2_eg_y - c2_gg_y;
  c2_we_b[1] = c2_ig_y + c2_kg_y;
  c2_we_b[2] = c2_lg_y - c2_mg_y;
  for (c2_i23 = 0; c2_i23 < 3; c2_i23++) {
    c2_dx2_dtheta[c2_i23] = c2_ig_a * c2_we_b[c2_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 78);
  for (c2_i24 = 0; c2_i24 < 9; c2_i24++) {
    c2_lf_a[c2_i24] = c2_T_e2E[c2_i24];
  }

  for (c2_i25 = 0; c2_i25 < 3; c2_i25++) {
    c2_we_b[c2_i25] = c2_force_mult_1[c2_i25];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i26 = 0; c2_i26 < 3; c2_i26++) {
    c2_mf_y[c2_i26] = 0.0;
    c2_i27 = 0;
    for (c2_i28 = 0; c2_i28 < 3; c2_i28++) {
      c2_mf_y[c2_i26] += c2_lf_a[c2_i27 + c2_i26] * c2_we_b[c2_i28];
      c2_i27 += 3;
    }
  }

  for (c2_i29 = 0; c2_i29 < 9; c2_i29++) {
    c2_lf_a[c2_i29] = c2_T_e2E[c2_i29];
  }

  for (c2_i30 = 0; c2_i30 < 3; c2_i30++) {
    c2_we_b[c2_i30] = c2_force_mult_2[c2_i30];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i31 = 0; c2_i31 < 3; c2_i31++) {
    c2_nf_y[c2_i31] = 0.0;
    c2_i32 = 0;
    for (c2_i33 = 0; c2_i33 < 3; c2_i33++) {
      c2_nf_y[c2_i31] += c2_lf_a[c2_i32 + c2_i31] * c2_we_b[c2_i33];
      c2_i32 += 3;
    }
  }

  for (c2_i34 = 0; c2_i34 < 3; c2_i34++) {
    c2_ng_y[c2_i34] = c2_mf_y[c2_i34];
  }

  for (c2_i35 = 0; c2_i35 < 3; c2_i35++) {
    c2_b_dx1_dtheta[c2_i35] = c2_dx1_dtheta[c2_i35];
  }

  for (c2_i36 = 0; c2_i36 < 3; c2_i36++) {
    c2_og_y[c2_i36] = c2_nf_y[c2_i36];
  }

  for (c2_i37 = 0; c2_i37 < 3; c2_i37++) {
    c2_b_dx2_dtheta[c2_i37] = c2_dx2_dtheta[c2_i37];
  }

  c2_upphi_theta = c2_dot(chartInstance, c2_ng_y, c2_b_dx1_dtheta) + c2_dot
    (chartInstance, c2_og_y, c2_b_dx2_dtheta);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 81);
  c2_e_A = c2_l;
  c2_lm_x = c2_e_A;
  c2_mm_x = c2_lm_x;
  c2_pg_y = c2_mm_x / 1.4142135623730951;
  c2_nm_x = c2_phi;
  c2_om_x = c2_nm_x;
  c2_om_x = muDoubleScalarCos(c2_om_x);
  c2_pm_x = c2_psi;
  c2_qm_x = c2_pm_x;
  c2_qm_x = muDoubleScalarSin(c2_qm_x);
  c2_jg_a = -c2_om_x;
  c2_fg_b = c2_qm_x;
  c2_qg_y = c2_jg_a * c2_fg_b;
  c2_rm_x = c2_theta;
  c2_sm_x = c2_rm_x;
  c2_sm_x = muDoubleScalarCos(c2_sm_x);
  c2_tm_x = c2_phi;
  c2_um_x = c2_tm_x;
  c2_um_x = muDoubleScalarSin(c2_um_x);
  c2_kg_a = c2_sm_x;
  c2_gg_b = c2_um_x;
  c2_rg_y = c2_kg_a * c2_gg_b;
  c2_vm_x = c2_psi;
  c2_wm_x = c2_vm_x;
  c2_wm_x = muDoubleScalarCos(c2_wm_x);
  c2_lg_a = c2_rg_y;
  c2_hg_b = c2_wm_x;
  c2_sg_y = c2_lg_a * c2_hg_b;
  c2_xm_x = c2_psi;
  c2_ym_x = c2_xm_x;
  c2_ym_x = muDoubleScalarCos(c2_ym_x);
  c2_an_x = c2_phi;
  c2_bn_x = c2_an_x;
  c2_bn_x = muDoubleScalarCos(c2_bn_x);
  c2_mg_a = c2_ym_x;
  c2_ig_b = c2_bn_x;
  c2_tg_y = c2_mg_a * c2_ig_b;
  c2_cn_x = c2_phi;
  c2_dn_x = c2_cn_x;
  c2_dn_x = muDoubleScalarSin(c2_dn_x);
  c2_en_x = c2_theta;
  c2_fn_x = c2_en_x;
  c2_fn_x = muDoubleScalarCos(c2_fn_x);
  c2_ng_a = c2_dn_x;
  c2_jg_b = c2_fn_x;
  c2_ug_y = c2_ng_a * c2_jg_b;
  c2_gn_x = c2_psi;
  c2_hn_x = c2_gn_x;
  c2_hn_x = muDoubleScalarSin(c2_hn_x);
  c2_og_a = c2_ug_y;
  c2_kg_b = c2_hn_x;
  c2_vg_y = c2_og_a * c2_kg_b;
  c2_in_x = c2_phi;
  c2_jn_x = c2_in_x;
  c2_jn_x = muDoubleScalarSin(c2_jn_x);
  c2_kn_x = c2_psi;
  c2_ln_x = c2_kn_x;
  c2_ln_x = muDoubleScalarSin(c2_ln_x);
  c2_pg_a = -c2_jn_x;
  c2_lg_b = c2_ln_x;
  c2_wg_y = c2_pg_a * c2_lg_b;
  c2_mn_x = c2_theta;
  c2_nn_x = c2_mn_x;
  c2_nn_x = muDoubleScalarCos(c2_nn_x);
  c2_on_x = c2_phi;
  c2_pn_x = c2_on_x;
  c2_pn_x = muDoubleScalarCos(c2_pn_x);
  c2_qg_a = c2_nn_x;
  c2_mg_b = c2_pn_x;
  c2_xg_y = c2_qg_a * c2_mg_b;
  c2_qn_x = c2_psi;
  c2_rn_x = c2_qn_x;
  c2_rn_x = muDoubleScalarCos(c2_rn_x);
  c2_rg_a = c2_xg_y;
  c2_ng_b = c2_rn_x;
  c2_yg_y = c2_rg_a * c2_ng_b;
  c2_sn_x = c2_phi;
  c2_tn_x = c2_sn_x;
  c2_tn_x = muDoubleScalarSin(c2_tn_x);
  c2_un_x = c2_psi;
  c2_vn_x = c2_un_x;
  c2_vn_x = muDoubleScalarCos(c2_vn_x);
  c2_sg_a = c2_tn_x;
  c2_og_b = c2_vn_x;
  c2_ah_y = c2_sg_a * c2_og_b;
  c2_wn_x = c2_phi;
  c2_xn_x = c2_wn_x;
  c2_xn_x = muDoubleScalarCos(c2_xn_x);
  c2_yn_x = c2_theta;
  c2_ao_x = c2_yn_x;
  c2_ao_x = muDoubleScalarCos(c2_ao_x);
  c2_tg_a = c2_xn_x;
  c2_pg_b = c2_ao_x;
  c2_bh_y = c2_tg_a * c2_pg_b;
  c2_bo_x = c2_psi;
  c2_co_x = c2_bo_x;
  c2_co_x = muDoubleScalarSin(c2_co_x);
  c2_ug_a = c2_bh_y;
  c2_qg_b = c2_co_x;
  c2_ch_y = c2_ug_a * c2_qg_b;
  c2_do_x = c2_theta;
  c2_eo_x = c2_do_x;
  c2_eo_x = muDoubleScalarSin(c2_eo_x);
  c2_fo_x = c2_psi;
  c2_go_x = c2_fo_x;
  c2_go_x = muDoubleScalarCos(c2_go_x);
  c2_vg_a = c2_eo_x;
  c2_rg_b = c2_go_x;
  c2_dh_y = c2_vg_a * c2_rg_b;
  c2_ho_x = c2_theta;
  c2_io_x = c2_ho_x;
  c2_io_x = muDoubleScalarSin(c2_io_x);
  c2_jo_x = c2_psi;
  c2_ko_x = c2_jo_x;
  c2_ko_x = muDoubleScalarSin(c2_ko_x);
  c2_wg_a = c2_io_x;
  c2_sg_b = c2_ko_x;
  c2_eh_y = c2_wg_a * c2_sg_b;
  c2_xg_a = c2_pg_y;
  c2_we_b[0] = ((c2_qg_y - c2_sg_y) - c2_tg_y) + c2_vg_y;
  c2_we_b[1] = ((c2_wg_y + c2_yg_y) - c2_ah_y) - c2_ch_y;
  c2_we_b[2] = c2_dh_y - c2_eh_y;
  for (c2_i38 = 0; c2_i38 < 3; c2_i38++) {
    c2_dx1_dpsi[c2_i38] = c2_xg_a * c2_we_b[c2_i38];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 86);
  c2_f_A = c2_l;
  c2_lo_x = c2_f_A;
  c2_mo_x = c2_lo_x;
  c2_fh_y = c2_mo_x / 1.4142135623730951;
  c2_no_x = c2_psi;
  c2_oo_x = c2_no_x;
  c2_oo_x = muDoubleScalarCos(c2_oo_x);
  c2_po_x = c2_phi;
  c2_qo_x = c2_po_x;
  c2_qo_x = muDoubleScalarCos(c2_qo_x);
  c2_yg_a = -c2_oo_x;
  c2_tg_b = c2_qo_x;
  c2_gh_y = c2_yg_a * c2_tg_b;
  c2_ro_x = c2_phi;
  c2_so_x = c2_ro_x;
  c2_so_x = muDoubleScalarSin(c2_so_x);
  c2_to_x = c2_theta;
  c2_uo_x = c2_to_x;
  c2_uo_x = muDoubleScalarCos(c2_uo_x);
  c2_ah_a = c2_so_x;
  c2_ug_b = c2_uo_x;
  c2_hh_y = c2_ah_a * c2_ug_b;
  c2_vo_x = c2_psi;
  c2_wo_x = c2_vo_x;
  c2_wo_x = muDoubleScalarSin(c2_wo_x);
  c2_bh_a = c2_hh_y;
  c2_vg_b = c2_wo_x;
  c2_ih_y = c2_bh_a * c2_vg_b;
  c2_xo_x = c2_phi;
  c2_yo_x = c2_xo_x;
  c2_yo_x = muDoubleScalarCos(c2_yo_x);
  c2_ap_x = c2_psi;
  c2_bp_x = c2_ap_x;
  c2_bp_x = muDoubleScalarSin(c2_bp_x);
  c2_ch_a = c2_yo_x;
  c2_wg_b = c2_bp_x;
  c2_jh_y = c2_ch_a * c2_wg_b;
  c2_cp_x = c2_theta;
  c2_dp_x = c2_cp_x;
  c2_dp_x = muDoubleScalarCos(c2_dp_x);
  c2_ep_x = c2_phi;
  c2_fp_x = c2_ep_x;
  c2_fp_x = muDoubleScalarSin(c2_fp_x);
  c2_dh_a = c2_dp_x;
  c2_xg_b = c2_fp_x;
  c2_kh_y = c2_dh_a * c2_xg_b;
  c2_gp_x = c2_psi;
  c2_hp_x = c2_gp_x;
  c2_hp_x = muDoubleScalarCos(c2_hp_x);
  c2_eh_a = c2_kh_y;
  c2_yg_b = c2_hp_x;
  c2_lh_y = c2_eh_a * c2_yg_b;
  c2_ip_x = c2_phi;
  c2_jp_x = c2_ip_x;
  c2_jp_x = muDoubleScalarSin(c2_jp_x);
  c2_kp_x = c2_psi;
  c2_lp_x = c2_kp_x;
  c2_lp_x = muDoubleScalarCos(c2_lp_x);
  c2_fh_a = -c2_jp_x;
  c2_ah_b = c2_lp_x;
  c2_mh_y = c2_fh_a * c2_ah_b;
  c2_mp_x = c2_phi;
  c2_np_x = c2_mp_x;
  c2_np_x = muDoubleScalarCos(c2_np_x);
  c2_op_x = c2_theta;
  c2_pp_x = c2_op_x;
  c2_pp_x = muDoubleScalarCos(c2_pp_x);
  c2_gh_a = c2_np_x;
  c2_bh_b = c2_pp_x;
  c2_nh_y = c2_gh_a * c2_bh_b;
  c2_qp_x = c2_psi;
  c2_rp_x = c2_qp_x;
  c2_rp_x = muDoubleScalarSin(c2_rp_x);
  c2_hh_a = c2_nh_y;
  c2_ch_b = c2_rp_x;
  c2_oh_y = c2_hh_a * c2_ch_b;
  c2_sp_x = c2_phi;
  c2_tp_x = c2_sp_x;
  c2_tp_x = muDoubleScalarSin(c2_tp_x);
  c2_up_x = c2_psi;
  c2_vp_x = c2_up_x;
  c2_vp_x = muDoubleScalarSin(c2_vp_x);
  c2_ih_a = c2_tp_x;
  c2_dh_b = c2_vp_x;
  c2_ph_y = c2_ih_a * c2_dh_b;
  c2_wp_x = c2_theta;
  c2_xp_x = c2_wp_x;
  c2_xp_x = muDoubleScalarCos(c2_xp_x);
  c2_yp_x = c2_phi;
  c2_aq_x = c2_yp_x;
  c2_aq_x = muDoubleScalarCos(c2_aq_x);
  c2_jh_a = c2_xp_x;
  c2_eh_b = c2_aq_x;
  c2_qh_y = c2_jh_a * c2_eh_b;
  c2_bq_x = c2_psi;
  c2_cq_x = c2_bq_x;
  c2_cq_x = muDoubleScalarCos(c2_cq_x);
  c2_kh_a = c2_qh_y;
  c2_fh_b = c2_cq_x;
  c2_rh_y = c2_kh_a * c2_fh_b;
  c2_dq_x = c2_theta;
  c2_eq_x = c2_dq_x;
  c2_eq_x = muDoubleScalarSin(c2_eq_x);
  c2_fq_x = c2_psi;
  c2_gq_x = c2_fq_x;
  c2_gq_x = muDoubleScalarSin(c2_gq_x);
  c2_lh_a = -c2_eq_x;
  c2_gh_b = c2_gq_x;
  c2_sh_y = c2_lh_a * c2_gh_b;
  c2_hq_x = c2_theta;
  c2_iq_x = c2_hq_x;
  c2_iq_x = muDoubleScalarSin(c2_iq_x);
  c2_jq_x = c2_psi;
  c2_kq_x = c2_jq_x;
  c2_kq_x = muDoubleScalarCos(c2_kq_x);
  c2_mh_a = c2_iq_x;
  c2_hh_b = c2_kq_x;
  c2_th_y = c2_mh_a * c2_hh_b;
  c2_nh_a = c2_fh_y;
  c2_we_b[0] = ((c2_gh_y + c2_ih_y) + c2_jh_y) + c2_lh_y;
  c2_we_b[1] = ((c2_mh_y - c2_oh_y) + c2_ph_y) - c2_rh_y;
  c2_we_b[2] = c2_sh_y - c2_th_y;
  for (c2_i39 = 0; c2_i39 < 3; c2_i39++) {
    c2_dx2_dpsi[c2_i39] = c2_nh_a * c2_we_b[c2_i39];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 91);
  for (c2_i40 = 0; c2_i40 < 9; c2_i40++) {
    c2_lf_a[c2_i40] = c2_T_e2E[c2_i40];
  }

  for (c2_i41 = 0; c2_i41 < 3; c2_i41++) {
    c2_we_b[c2_i41] = c2_force_mult_1[c2_i41];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i42 = 0; c2_i42 < 3; c2_i42++) {
    c2_mf_y[c2_i42] = 0.0;
    c2_i43 = 0;
    for (c2_i44 = 0; c2_i44 < 3; c2_i44++) {
      c2_mf_y[c2_i42] += c2_lf_a[c2_i43 + c2_i42] * c2_we_b[c2_i44];
      c2_i43 += 3;
    }
  }

  for (c2_i45 = 0; c2_i45 < 9; c2_i45++) {
    c2_lf_a[c2_i45] = c2_T_e2E[c2_i45];
  }

  for (c2_i46 = 0; c2_i46 < 3; c2_i46++) {
    c2_we_b[c2_i46] = c2_force_mult_2[c2_i46];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i47 = 0; c2_i47 < 3; c2_i47++) {
    c2_nf_y[c2_i47] = 0.0;
    c2_i48 = 0;
    for (c2_i49 = 0; c2_i49 < 3; c2_i49++) {
      c2_nf_y[c2_i47] += c2_lf_a[c2_i48 + c2_i47] * c2_we_b[c2_i49];
      c2_i48 += 3;
    }
  }

  for (c2_i50 = 0; c2_i50 < 3; c2_i50++) {
    c2_uh_y[c2_i50] = c2_mf_y[c2_i50];
  }

  for (c2_i51 = 0; c2_i51 < 3; c2_i51++) {
    c2_b_dx1_dpsi[c2_i51] = c2_dx1_dpsi[c2_i51];
  }

  for (c2_i52 = 0; c2_i52 < 3; c2_i52++) {
    c2_vh_y[c2_i52] = c2_nf_y[c2_i52];
  }

  for (c2_i53 = 0; c2_i53 < 3; c2_i53++) {
    c2_b_dx2_dpsi[c2_i53] = c2_dx2_dpsi[c2_i53];
  }

  c2_upphi_psi = c2_dot(chartInstance, c2_uh_y, c2_b_dx1_dpsi) + c2_dot
    (chartInstance, c2_vh_y, c2_b_dx2_dpsi);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 93);
  c2_ih_b = c2_Jxx;
  c2_wh_y = 2.0 * c2_ih_b;
  c2_lq_x = c2_psi;
  c2_mq_x = c2_lq_x;
  c2_mq_x = muDoubleScalarSin(c2_mq_x);
  c2_oh_a = c2_wh_y;
  c2_jh_b = c2_mpower(chartInstance, c2_mq_x);
  c2_xh_y = c2_oh_a * c2_jh_b;
  c2_nq_x = c2_theta;
  c2_oq_x = c2_nq_x;
  c2_oq_x = muDoubleScalarSin(c2_oq_x);
  c2_ph_a = c2_xh_y;
  c2_kh_b = c2_oq_x;
  c2_yh_y = c2_ph_a * c2_kh_b;
  c2_pq_x = c2_theta;
  c2_qq_x = c2_pq_x;
  c2_qq_x = muDoubleScalarCos(c2_qq_x);
  c2_qh_a = c2_yh_y;
  c2_lh_b = c2_qq_x;
  c2_ai_y = c2_qh_a * c2_lh_b;
  c2_rh_a = c2_ai_y;
  c2_mh_b = c2_theta_dot;
  c2_bi_y = c2_rh_a * c2_mh_b;
  c2_sh_a = c2_bi_y;
  c2_nh_b = c2_phi_dot;
  c2_ci_y = c2_sh_a * c2_nh_b;
  c2_oh_b = c2_Jxx;
  c2_di_y = 2.0 * c2_oh_b;
  c2_rq_x = c2_psi;
  c2_sq_x = c2_rq_x;
  c2_sq_x = muDoubleScalarSin(c2_sq_x);
  c2_th_a = c2_di_y;
  c2_ph_b = c2_sq_x;
  c2_ei_y = c2_th_a * c2_ph_b;
  c2_tq_x = c2_psi;
  c2_uq_x = c2_tq_x;
  c2_uq_x = muDoubleScalarCos(c2_uq_x);
  c2_uh_a = c2_ei_y;
  c2_qh_b = c2_uq_x;
  c2_fi_y = c2_uh_a * c2_qh_b;
  c2_vq_x = c2_theta;
  c2_wq_x = c2_vq_x;
  c2_wq_x = muDoubleScalarSin(c2_wq_x);
  c2_vh_a = c2_fi_y;
  c2_rh_b = c2_mpower(chartInstance, c2_wq_x);
  c2_gi_y = c2_vh_a * c2_rh_b;
  c2_wh_a = c2_gi_y;
  c2_sh_b = c2_psi_dot;
  c2_hi_y = c2_wh_a * c2_sh_b;
  c2_xh_a = c2_hi_y;
  c2_th_b = c2_phi_dot;
  c2_ii_y = c2_xh_a * c2_th_b;
  c2_xq_x = c2_psi;
  c2_yq_x = c2_xq_x;
  c2_yq_x = muDoubleScalarSin(c2_yq_x);
  c2_yh_a = c2_Jxx;
  c2_uh_b = c2_mpower(chartInstance, c2_yq_x);
  c2_ji_y = c2_yh_a * c2_uh_b;
  c2_ar_x = c2_theta;
  c2_br_x = c2_ar_x;
  c2_br_x = muDoubleScalarSin(c2_br_x);
  c2_ai_a = c2_ji_y;
  c2_vh_b = c2_br_x;
  c2_ki_y = c2_ai_a * c2_vh_b;
  c2_bi_a = c2_ki_y;
  c2_wh_b = c2_psi_dot;
  c2_li_y = c2_bi_a * c2_wh_b;
  c2_ci_a = c2_li_y;
  c2_xh_b = c2_theta_dot;
  c2_mi_y = c2_ci_a * c2_xh_b;
  c2_cr_x = c2_psi;
  c2_dr_x = c2_cr_x;
  c2_dr_x = muDoubleScalarSin(c2_dr_x);
  c2_di_a = c2_Jxx;
  c2_yh_b = c2_dr_x;
  c2_ni_y = c2_di_a * c2_yh_b;
  c2_er_x = c2_theta;
  c2_fr_x = c2_er_x;
  c2_fr_x = muDoubleScalarCos(c2_fr_x);
  c2_ei_a = c2_ni_y;
  c2_ai_b = c2_fr_x;
  c2_oi_y = c2_ei_a * c2_ai_b;
  c2_gr_x = c2_psi;
  c2_hr_x = c2_gr_x;
  c2_hr_x = muDoubleScalarCos(c2_hr_x);
  c2_fi_a = c2_oi_y;
  c2_bi_b = c2_hr_x;
  c2_pi_y = c2_fi_a * c2_bi_b;
  c2_gi_a = c2_pi_y;
  c2_ci_b = c2_mpower(chartInstance, c2_theta_dot);
  c2_qi_y = c2_gi_a * c2_ci_b;
  c2_ir_x = c2_psi;
  c2_jr_x = c2_ir_x;
  c2_jr_x = muDoubleScalarCos(c2_jr_x);
  c2_hi_a = c2_Jxx;
  c2_di_b = c2_mpower(chartInstance, c2_jr_x);
  c2_ri_y = c2_hi_a * c2_di_b;
  c2_kr_x = c2_theta;
  c2_lr_x = c2_kr_x;
  c2_lr_x = muDoubleScalarSin(c2_lr_x);
  c2_ii_a = c2_ri_y;
  c2_ei_b = c2_lr_x;
  c2_si_y = c2_ii_a * c2_ei_b;
  c2_ji_a = c2_si_y;
  c2_fi_b = c2_psi_dot;
  c2_ti_y = c2_ji_a * c2_fi_b;
  c2_ki_a = c2_ti_y;
  c2_gi_b = c2_theta_dot;
  c2_ui_y = c2_ki_a * c2_gi_b;
  c2_hi_b = c2_Jyy;
  c2_vi_y = 2.0 * c2_hi_b;
  c2_mr_x = c2_theta;
  c2_nr_x = c2_mr_x;
  c2_nr_x = muDoubleScalarSin(c2_nr_x);
  c2_li_a = c2_vi_y;
  c2_ii_b = c2_nr_x;
  c2_wi_y = c2_li_a * c2_ii_b;
  c2_or_x = c2_theta;
  c2_pr_x = c2_or_x;
  c2_pr_x = muDoubleScalarCos(c2_pr_x);
  c2_mi_a = c2_wi_y;
  c2_ji_b = c2_pr_x;
  c2_xi_y = c2_mi_a * c2_ji_b;
  c2_qr_x = c2_psi;
  c2_rr_x = c2_qr_x;
  c2_rr_x = muDoubleScalarCos(c2_rr_x);
  c2_ni_a = c2_xi_y;
  c2_ki_b = c2_mpower(chartInstance, c2_rr_x);
  c2_yi_y = c2_ni_a * c2_ki_b;
  c2_oi_a = c2_yi_y;
  c2_li_b = c2_theta_dot;
  c2_aj_y = c2_oi_a * c2_li_b;
  c2_pi_a = c2_aj_y;
  c2_mi_b = c2_phi_dot;
  c2_bj_y = c2_pi_a * c2_mi_b;
  c2_ni_b = c2_Jyy;
  c2_cj_y = 2.0 * c2_ni_b;
  c2_sr_x = c2_psi;
  c2_tr_x = c2_sr_x;
  c2_tr_x = muDoubleScalarCos(c2_tr_x);
  c2_qi_a = c2_cj_y;
  c2_oi_b = c2_tr_x;
  c2_dj_y = c2_qi_a * c2_oi_b;
  c2_ur_x = c2_psi;
  c2_vr_x = c2_ur_x;
  c2_vr_x = muDoubleScalarSin(c2_vr_x);
  c2_ri_a = c2_dj_y;
  c2_pi_b = c2_vr_x;
  c2_ej_y = c2_ri_a * c2_pi_b;
  c2_wr_x = c2_theta;
  c2_xr_x = c2_wr_x;
  c2_xr_x = muDoubleScalarSin(c2_xr_x);
  c2_si_a = c2_ej_y;
  c2_qi_b = c2_mpower(chartInstance, c2_xr_x);
  c2_fj_y = c2_si_a * c2_qi_b;
  c2_ti_a = c2_fj_y;
  c2_ri_b = c2_psi_dot;
  c2_gj_y = c2_ti_a * c2_ri_b;
  c2_ui_a = c2_gj_y;
  c2_si_b = c2_phi_dot;
  c2_hj_y = c2_ui_a * c2_si_b;
  c2_yr_x = c2_psi;
  c2_as_x = c2_yr_x;
  c2_as_x = muDoubleScalarCos(c2_as_x);
  c2_vi_a = c2_Jyy;
  c2_ti_b = c2_mpower(chartInstance, c2_as_x);
  c2_ij_y = c2_vi_a * c2_ti_b;
  c2_bs_x = c2_theta;
  c2_cs_x = c2_bs_x;
  c2_cs_x = muDoubleScalarSin(c2_cs_x);
  c2_wi_a = c2_ij_y;
  c2_ui_b = c2_cs_x;
  c2_jj_y = c2_wi_a * c2_ui_b;
  c2_xi_a = c2_jj_y;
  c2_vi_b = c2_psi_dot;
  c2_kj_y = c2_xi_a * c2_vi_b;
  c2_yi_a = c2_kj_y;
  c2_wi_b = c2_theta_dot;
  c2_lj_y = c2_yi_a * c2_wi_b;
  c2_ds_x = c2_psi;
  c2_es_x = c2_ds_x;
  c2_es_x = muDoubleScalarCos(c2_es_x);
  c2_aj_a = c2_Jyy;
  c2_xi_b = c2_es_x;
  c2_mj_y = c2_aj_a * c2_xi_b;
  c2_fs_x = c2_theta;
  c2_gs_x = c2_fs_x;
  c2_gs_x = muDoubleScalarCos(c2_gs_x);
  c2_bj_a = c2_mj_y;
  c2_yi_b = c2_gs_x;
  c2_nj_y = c2_bj_a * c2_yi_b;
  c2_hs_x = c2_psi;
  c2_is_x = c2_hs_x;
  c2_is_x = muDoubleScalarSin(c2_is_x);
  c2_cj_a = c2_nj_y;
  c2_aj_b = c2_is_x;
  c2_oj_y = c2_cj_a * c2_aj_b;
  c2_dj_a = c2_oj_y;
  c2_bj_b = c2_mpower(chartInstance, c2_theta_dot);
  c2_pj_y = c2_dj_a * c2_bj_b;
  c2_js_x = c2_psi;
  c2_ks_x = c2_js_x;
  c2_ks_x = muDoubleScalarSin(c2_ks_x);
  c2_ej_a = c2_Jyy;
  c2_cj_b = c2_mpower(chartInstance, c2_ks_x);
  c2_qj_y = c2_ej_a * c2_cj_b;
  c2_ls_x = c2_theta;
  c2_ms_x = c2_ls_x;
  c2_ms_x = muDoubleScalarSin(c2_ms_x);
  c2_fj_a = c2_qj_y;
  c2_dj_b = c2_ms_x;
  c2_rj_y = c2_fj_a * c2_dj_b;
  c2_gj_a = c2_rj_y;
  c2_ej_b = c2_psi_dot;
  c2_sj_y = c2_gj_a * c2_ej_b;
  c2_hj_a = c2_sj_y;
  c2_fj_b = c2_theta_dot;
  c2_tj_y = c2_hj_a * c2_fj_b;
  c2_gj_b = c2_Jzz;
  c2_uj_y = 2.0 * c2_gj_b;
  c2_ns_x = c2_theta;
  c2_os_x = c2_ns_x;
  c2_os_x = muDoubleScalarCos(c2_os_x);
  c2_ij_a = c2_uj_y;
  c2_hj_b = c2_os_x;
  c2_vj_y = c2_ij_a * c2_hj_b;
  c2_ps_x = c2_theta;
  c2_qs_x = c2_ps_x;
  c2_qs_x = muDoubleScalarSin(c2_qs_x);
  c2_jj_a = c2_vj_y;
  c2_ij_b = c2_qs_x;
  c2_wj_y = c2_jj_a * c2_ij_b;
  c2_kj_a = c2_wj_y;
  c2_jj_b = c2_theta_dot;
  c2_xj_y = c2_kj_a * c2_jj_b;
  c2_lj_a = c2_xj_y;
  c2_kj_b = c2_phi_dot;
  c2_yj_y = c2_lj_a * c2_kj_b;
  c2_rs_x = c2_theta;
  c2_ss_x = c2_rs_x;
  c2_ss_x = muDoubleScalarSin(c2_ss_x);
  c2_mj_a = c2_Jzz;
  c2_lj_b = c2_ss_x;
  c2_ak_y = c2_mj_a * c2_lj_b;
  c2_nj_a = c2_ak_y;
  c2_mj_b = c2_theta_dot;
  c2_bk_y = c2_nj_a * c2_mj_b;
  c2_oj_a = c2_bk_y;
  c2_nj_b = c2_psi_dot;
  c2_ck_y = c2_oj_a * c2_nj_b;
  c2_delta1 = (((((((((((c2_upphi_phi - c2_ci_y) - c2_ii_y) + c2_mi_y) - c2_qi_y)
                     - c2_ui_y) - c2_bj_y) + c2_hj_y) + c2_lj_y) + c2_pj_y) -
                c2_tj_y) + c2_yj_y) + c2_ck_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 100);
  c2_ts_x = c2_psi;
  c2_us_x = c2_ts_x;
  c2_us_x = muDoubleScalarSin(c2_us_x);
  c2_pj_a = c2_Jxx;
  c2_oj_b = c2_mpower(chartInstance, c2_us_x);
  c2_dk_y = c2_pj_a * c2_oj_b;
  c2_vs_x = c2_theta;
  c2_ws_x = c2_vs_x;
  c2_ws_x = muDoubleScalarSin(c2_ws_x);
  c2_qj_a = c2_dk_y;
  c2_pj_b = c2_ws_x;
  c2_ek_y = c2_qj_a * c2_pj_b;
  c2_rj_a = c2_ek_y;
  c2_qj_b = c2_psi_dot;
  c2_fk_y = c2_rj_a * c2_qj_b;
  c2_sj_a = c2_fk_y;
  c2_rj_b = c2_phi_dot;
  c2_gk_y = c2_sj_a * c2_rj_b;
  c2_xs_x = c2_psi;
  c2_ys_x = c2_xs_x;
  c2_ys_x = muDoubleScalarCos(c2_ys_x);
  c2_tj_a = c2_Jxx;
  c2_sj_b = c2_mpower(chartInstance, c2_ys_x);
  c2_hk_y = c2_tj_a * c2_sj_b;
  c2_at_x = c2_theta;
  c2_bt_x = c2_at_x;
  c2_bt_x = muDoubleScalarSin(c2_bt_x);
  c2_uj_a = c2_hk_y;
  c2_tj_b = c2_bt_x;
  c2_ik_y = c2_uj_a * c2_tj_b;
  c2_vj_a = c2_ik_y;
  c2_uj_b = c2_psi_dot;
  c2_jk_y = c2_vj_a * c2_uj_b;
  c2_wj_a = c2_jk_y;
  c2_vj_b = c2_phi_dot;
  c2_kk_y = c2_wj_a * c2_vj_b;
  c2_wj_b = c2_Jxx;
  c2_lk_y = 2.0 * c2_wj_b;
  c2_ct_x = c2_psi;
  c2_dt_x = c2_ct_x;
  c2_dt_x = muDoubleScalarCos(c2_dt_x);
  c2_xj_a = c2_lk_y;
  c2_xj_b = c2_dt_x;
  c2_mk_y = c2_xj_a * c2_xj_b;
  c2_et_x = c2_psi;
  c2_ft_x = c2_et_x;
  c2_ft_x = muDoubleScalarSin(c2_ft_x);
  c2_yj_a = c2_mk_y;
  c2_yj_b = c2_ft_x;
  c2_nk_y = c2_yj_a * c2_yj_b;
  c2_ak_a = c2_nk_y;
  c2_ak_b = c2_psi_dot;
  c2_ok_y = c2_ak_a * c2_ak_b;
  c2_bk_a = c2_ok_y;
  c2_bk_b = c2_theta_dot;
  c2_pk_y = c2_bk_a * c2_bk_b;
  c2_gt_x = c2_psi;
  c2_ht_x = c2_gt_x;
  c2_ht_x = muDoubleScalarSin(c2_ht_x);
  c2_ck_a = c2_Jxx;
  c2_ck_b = c2_mpower(chartInstance, c2_ht_x);
  c2_qk_y = c2_ck_a * c2_ck_b;
  c2_it_x = c2_theta;
  c2_jt_x = c2_it_x;
  c2_jt_x = muDoubleScalarSin(c2_jt_x);
  c2_dk_a = c2_qk_y;
  c2_dk_b = c2_jt_x;
  c2_rk_y = c2_dk_a * c2_dk_b;
  c2_kt_x = c2_theta;
  c2_lt_x = c2_kt_x;
  c2_lt_x = muDoubleScalarCos(c2_lt_x);
  c2_ek_a = c2_rk_y;
  c2_ek_b = c2_lt_x;
  c2_sk_y = c2_ek_a * c2_ek_b;
  c2_fk_a = c2_sk_y;
  c2_fk_b = c2_mpower(chartInstance, c2_phi_dot);
  c2_tk_y = c2_fk_a * c2_fk_b;
  c2_mt_x = c2_psi;
  c2_nt_x = c2_mt_x;
  c2_nt_x = muDoubleScalarCos(c2_nt_x);
  c2_gk_a = c2_Jyy;
  c2_gk_b = c2_mpower(chartInstance, c2_nt_x);
  c2_uk_y = c2_gk_a * c2_gk_b;
  c2_ot_x = c2_theta;
  c2_pt_x = c2_ot_x;
  c2_pt_x = muDoubleScalarSin(c2_pt_x);
  c2_hk_a = c2_uk_y;
  c2_hk_b = c2_pt_x;
  c2_vk_y = c2_hk_a * c2_hk_b;
  c2_ik_a = c2_vk_y;
  c2_ik_b = c2_psi_dot;
  c2_wk_y = c2_ik_a * c2_ik_b;
  c2_jk_a = c2_wk_y;
  c2_jk_b = c2_phi_dot;
  c2_xk_y = c2_jk_a * c2_jk_b;
  c2_qt_x = c2_psi;
  c2_rt_x = c2_qt_x;
  c2_rt_x = muDoubleScalarSin(c2_rt_x);
  c2_kk_a = c2_Jyy;
  c2_kk_b = c2_mpower(chartInstance, c2_rt_x);
  c2_yk_y = c2_kk_a * c2_kk_b;
  c2_st_x = c2_theta;
  c2_tt_x = c2_st_x;
  c2_tt_x = muDoubleScalarSin(c2_tt_x);
  c2_lk_a = c2_yk_y;
  c2_lk_b = c2_tt_x;
  c2_al_y = c2_lk_a * c2_lk_b;
  c2_mk_a = c2_al_y;
  c2_mk_b = c2_psi_dot;
  c2_bl_y = c2_mk_a * c2_mk_b;
  c2_nk_a = c2_bl_y;
  c2_nk_b = c2_phi_dot;
  c2_cl_y = c2_nk_a * c2_nk_b;
  c2_ok_b = c2_Jyy;
  c2_dl_y = 2.0 * c2_ok_b;
  c2_ut_x = c2_psi;
  c2_vt_x = c2_ut_x;
  c2_vt_x = muDoubleScalarSin(c2_vt_x);
  c2_ok_a = c2_dl_y;
  c2_pk_b = c2_vt_x;
  c2_el_y = c2_ok_a * c2_pk_b;
  c2_wt_x = c2_psi;
  c2_xt_x = c2_wt_x;
  c2_xt_x = muDoubleScalarCos(c2_xt_x);
  c2_pk_a = c2_el_y;
  c2_qk_b = c2_xt_x;
  c2_fl_y = c2_pk_a * c2_qk_b;
  c2_qk_a = c2_fl_y;
  c2_rk_b = c2_psi_dot;
  c2_gl_y = c2_qk_a * c2_rk_b;
  c2_rk_a = c2_gl_y;
  c2_sk_b = c2_theta_dot;
  c2_hl_y = c2_rk_a * c2_sk_b;
  c2_yt_x = c2_psi;
  c2_au_x = c2_yt_x;
  c2_au_x = muDoubleScalarCos(c2_au_x);
  c2_sk_a = c2_Jyy;
  c2_tk_b = c2_mpower(chartInstance, c2_au_x);
  c2_il_y = c2_sk_a * c2_tk_b;
  c2_bu_x = c2_theta;
  c2_cu_x = c2_bu_x;
  c2_cu_x = muDoubleScalarSin(c2_cu_x);
  c2_tk_a = c2_il_y;
  c2_uk_b = c2_cu_x;
  c2_jl_y = c2_tk_a * c2_uk_b;
  c2_du_x = c2_theta;
  c2_eu_x = c2_du_x;
  c2_eu_x = muDoubleScalarCos(c2_eu_x);
  c2_uk_a = c2_jl_y;
  c2_vk_b = c2_eu_x;
  c2_kl_y = c2_uk_a * c2_vk_b;
  c2_vk_a = c2_kl_y;
  c2_wk_b = c2_mpower(chartInstance, c2_phi_dot);
  c2_ll_y = c2_vk_a * c2_wk_b;
  c2_fu_x = c2_theta;
  c2_gu_x = c2_fu_x;
  c2_gu_x = muDoubleScalarSin(c2_gu_x);
  c2_wk_a = c2_Jzz;
  c2_xk_b = c2_gu_x;
  c2_ml_y = c2_wk_a * c2_xk_b;
  c2_hu_x = c2_theta;
  c2_iu_x = c2_hu_x;
  c2_iu_x = muDoubleScalarCos(c2_iu_x);
  c2_xk_a = c2_ml_y;
  c2_yk_b = c2_iu_x;
  c2_nl_y = c2_xk_a * c2_yk_b;
  c2_yk_a = c2_nl_y;
  c2_al_b = c2_mpower(chartInstance, c2_phi_dot);
  c2_ol_y = c2_yk_a * c2_al_b;
  c2_ju_x = c2_theta;
  c2_ku_x = c2_ju_x;
  c2_ku_x = muDoubleScalarSin(c2_ku_x);
  c2_al_a = c2_Jzz;
  c2_bl_b = c2_ku_x;
  c2_pl_y = c2_al_a * c2_bl_b;
  c2_bl_a = c2_pl_y;
  c2_cl_b = c2_phi_dot;
  c2_ql_y = c2_bl_a * c2_cl_b;
  c2_cl_a = c2_ql_y;
  c2_dl_b = c2_psi_dot;
  c2_rl_y = c2_cl_a * c2_dl_b;
  c2_delta2 = (((((((((c2_upphi_theta + c2_gk_y) - c2_kk_y) + c2_pk_y) + c2_tk_y)
                   + c2_xk_y) - c2_cl_y) - c2_hl_y) + c2_ll_y) - c2_ol_y) -
    c2_rl_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 107);
  c2_lu_x = c2_theta;
  c2_mu_x = c2_lu_x;
  c2_mu_x = muDoubleScalarSin(c2_mu_x);
  c2_dl_a = c2_Jxx;
  c2_el_b = c2_mpower(chartInstance, c2_mu_x);
  c2_sl_y = c2_dl_a * c2_el_b;
  c2_nu_x = c2_psi;
  c2_ou_x = c2_nu_x;
  c2_ou_x = muDoubleScalarSin(c2_ou_x);
  c2_el_a = c2_sl_y;
  c2_fl_b = c2_ou_x;
  c2_tl_y = c2_el_a * c2_fl_b;
  c2_pu_x = c2_psi;
  c2_qu_x = c2_pu_x;
  c2_qu_x = muDoubleScalarCos(c2_qu_x);
  c2_fl_a = c2_tl_y;
  c2_gl_b = c2_qu_x;
  c2_ul_y = c2_fl_a * c2_gl_b;
  c2_gl_a = c2_ul_y;
  c2_hl_b = c2_mpower(chartInstance, c2_phi_dot);
  c2_vl_y = c2_gl_a * c2_hl_b;
  c2_ru_x = c2_theta;
  c2_su_x = c2_ru_x;
  c2_su_x = muDoubleScalarSin(c2_su_x);
  c2_hl_a = c2_Jxx;
  c2_il_b = c2_su_x;
  c2_wl_y = c2_hl_a * c2_il_b;
  c2_tu_x = c2_psi;
  c2_uu_x = c2_tu_x;
  c2_uu_x = muDoubleScalarSin(c2_uu_x);
  c2_il_a = c2_wl_y;
  c2_jl_b = c2_mpower(chartInstance, c2_uu_x);
  c2_xl_y = c2_il_a * c2_jl_b;
  c2_jl_a = c2_xl_y;
  c2_kl_b = c2_phi_dot;
  c2_yl_y = c2_jl_a * c2_kl_b;
  c2_kl_a = c2_yl_y;
  c2_ll_b = c2_theta_dot;
  c2_am_y = c2_kl_a * c2_ll_b;
  c2_vu_x = c2_theta;
  c2_wu_x = c2_vu_x;
  c2_wu_x = muDoubleScalarSin(c2_wu_x);
  c2_ll_a = c2_Jxx;
  c2_ml_b = c2_wu_x;
  c2_bm_y = c2_ll_a * c2_ml_b;
  c2_xu_x = c2_psi;
  c2_yu_x = c2_xu_x;
  c2_yu_x = muDoubleScalarCos(c2_yu_x);
  c2_ml_a = c2_bm_y;
  c2_nl_b = c2_mpower(chartInstance, c2_yu_x);
  c2_cm_y = c2_ml_a * c2_nl_b;
  c2_nl_a = c2_cm_y;
  c2_ol_b = c2_phi_dot;
  c2_dm_y = c2_nl_a * c2_ol_b;
  c2_ol_a = c2_dm_y;
  c2_pl_b = c2_theta_dot;
  c2_em_y = c2_ol_a * c2_pl_b;
  c2_av_x = c2_psi;
  c2_bv_x = c2_av_x;
  c2_bv_x = muDoubleScalarCos(c2_bv_x);
  c2_pl_a = c2_Jxx;
  c2_ql_b = c2_bv_x;
  c2_fm_y = c2_pl_a * c2_ql_b;
  c2_cv_x = c2_psi;
  c2_dv_x = c2_cv_x;
  c2_dv_x = muDoubleScalarSin(c2_dv_x);
  c2_ql_a = c2_fm_y;
  c2_rl_b = c2_dv_x;
  c2_gm_y = c2_ql_a * c2_rl_b;
  c2_rl_a = c2_gm_y;
  c2_sl_b = c2_mpower(chartInstance, c2_theta_dot);
  c2_hm_y = c2_rl_a * c2_sl_b;
  c2_ev_x = c2_theta;
  c2_fv_x = c2_ev_x;
  c2_fv_x = muDoubleScalarSin(c2_fv_x);
  c2_sl_a = c2_Jyy;
  c2_tl_b = c2_mpower(chartInstance, c2_fv_x);
  c2_im_y = c2_sl_a * c2_tl_b;
  c2_gv_x = c2_psi;
  c2_hv_x = c2_gv_x;
  c2_hv_x = muDoubleScalarCos(c2_hv_x);
  c2_tl_a = c2_im_y;
  c2_ul_b = c2_hv_x;
  c2_jm_y = c2_tl_a * c2_ul_b;
  c2_iv_x = c2_psi;
  c2_jv_x = c2_iv_x;
  c2_jv_x = muDoubleScalarSin(c2_jv_x);
  c2_ul_a = c2_jm_y;
  c2_vl_b = c2_jv_x;
  c2_km_y = c2_ul_a * c2_vl_b;
  c2_vl_a = c2_km_y;
  c2_wl_b = c2_mpower(chartInstance, c2_phi_dot);
  c2_lm_y = c2_vl_a * c2_wl_b;
  c2_kv_x = c2_theta;
  c2_lv_x = c2_kv_x;
  c2_lv_x = muDoubleScalarSin(c2_lv_x);
  c2_wl_a = c2_Jyy;
  c2_xl_b = c2_lv_x;
  c2_mm_y = c2_wl_a * c2_xl_b;
  c2_mv_x = c2_psi;
  c2_nv_x = c2_mv_x;
  c2_nv_x = muDoubleScalarCos(c2_nv_x);
  c2_xl_a = c2_mm_y;
  c2_yl_b = c2_mpower(chartInstance, c2_nv_x);
  c2_nm_y = c2_xl_a * c2_yl_b;
  c2_yl_a = c2_nm_y;
  c2_am_b = c2_phi_dot;
  c2_om_y = c2_yl_a * c2_am_b;
  c2_am_a = c2_om_y;
  c2_bm_b = c2_theta_dot;
  c2_pm_y = c2_am_a * c2_bm_b;
  c2_ov_x = c2_theta;
  c2_pv_x = c2_ov_x;
  c2_pv_x = muDoubleScalarSin(c2_pv_x);
  c2_bm_a = c2_Jyy;
  c2_cm_b = c2_pv_x;
  c2_qm_y = c2_bm_a * c2_cm_b;
  c2_qv_x = c2_psi;
  c2_rv_x = c2_qv_x;
  c2_rv_x = muDoubleScalarSin(c2_rv_x);
  c2_cm_a = c2_qm_y;
  c2_dm_b = c2_mpower(chartInstance, c2_rv_x);
  c2_rm_y = c2_cm_a * c2_dm_b;
  c2_dm_a = c2_rm_y;
  c2_em_b = c2_phi_dot;
  c2_sm_y = c2_dm_a * c2_em_b;
  c2_em_a = c2_sm_y;
  c2_fm_b = c2_theta_dot;
  c2_tm_y = c2_em_a * c2_fm_b;
  c2_sv_x = c2_psi;
  c2_tv_x = c2_sv_x;
  c2_tv_x = muDoubleScalarSin(c2_tv_x);
  c2_fm_a = c2_Jyy;
  c2_gm_b = c2_tv_x;
  c2_um_y = c2_fm_a * c2_gm_b;
  c2_uv_x = c2_psi;
  c2_vv_x = c2_uv_x;
  c2_vv_x = muDoubleScalarCos(c2_vv_x);
  c2_gm_a = c2_um_y;
  c2_hm_b = c2_vv_x;
  c2_vm_y = c2_gm_a * c2_hm_b;
  c2_hm_a = c2_vm_y;
  c2_im_b = c2_mpower(chartInstance, c2_theta_dot);
  c2_wm_y = c2_hm_a * c2_im_b;
  c2_wv_x = c2_theta;
  c2_xv_x = c2_wv_x;
  c2_xv_x = muDoubleScalarSin(c2_xv_x);
  c2_im_a = c2_Jzz;
  c2_jm_b = c2_xv_x;
  c2_xm_y = c2_im_a * c2_jm_b;
  c2_jm_a = c2_xm_y;
  c2_km_b = c2_theta_dot;
  c2_ym_y = c2_jm_a * c2_km_b;
  c2_km_a = c2_ym_y;
  c2_lm_b = c2_phi_dot;
  c2_an_y = c2_km_a * c2_lm_b;
  c2_delta3 = ((((((((c2_upphi_psi + c2_vl_y) - c2_am_y) + c2_em_y) - c2_hm_y) -
                  c2_lm_y) - c2_pm_y) + c2_tm_y) + c2_wm_y) + c2_an_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 114);
  for (c2_i54 = 0; c2_i54 < 9; c2_i54++) {
    c2_b_multipliers[c2_i54] = c2_multipliers[c2_i54];
  }

  c2_b_delta1[0] = c2_delta1;
  c2_b_delta1[1] = c2_delta2;
  c2_b_delta1[2] = c2_delta3;
  c2_mldivide(chartInstance, c2_b_multipliers, c2_b_delta1, c2_dv0);
  for (c2_i55 = 0; c2_i55 < 3; c2_i55++) {
    c2_ddot[c2_i55] = c2_dv0[c2_i55];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 117);
  c2_g_A = c2_upphi_x;
  c2_B = c2_m;
  c2_yv_x = c2_g_A;
  c2_bn_y = c2_B;
  c2_aw_x = c2_yv_x;
  c2_cn_y = c2_bn_y;
  c2_x_ddot = c2_aw_x / c2_cn_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 118);
  c2_h_A = c2_upphi_y;
  c2_b_B = c2_m;
  c2_bw_x = c2_h_A;
  c2_dn_y = c2_b_B;
  c2_cw_x = c2_bw_x;
  c2_en_y = c2_dn_y;
  c2_y_ddot = c2_cw_x / c2_en_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 119);
  c2_i_A = c2_upphi_z;
  c2_c_B = c2_m;
  c2_dw_x = c2_i_A;
  c2_fn_y = c2_c_B;
  c2_ew_x = c2_dw_x;
  c2_gn_y = c2_fn_y;
  c2_hn_y = c2_ew_x / c2_gn_y;
  c2_z_ddot = c2_hn_y - c2_g;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 120);
  c2_phi_ddot = c2_ddot[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 121);
  c2_theta_ddot = c2_ddot[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 122);
  c2_psi_ddot = c2_ddot[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -122);
  _SFD_SYMBOL_SCOPE_POP();
  *c2_b_x_ddot = c2_x_ddot;
  *c2_b_y_ddot = c2_y_ddot;
  *c2_b_z_ddot = c2_z_ddot;
  *c2_b_phi_ddot = c2_phi_ddot;
  *c2_b_theta_ddot = c2_theta_ddot;
  *c2_b_psi_ddot = c2_psi_ddot;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c2_sfEvent);
}

static void initSimStructsc2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
}

static void registerMessagesc2_quad_dynamics_sim_e
  (SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber)
{
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static real_T c2_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_psi_ddot, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_psi_ddot),
    &c2_thisId);
  sf_mex_destroy(&c2_psi_ddot);
  return c2_y;
}

static real_T c2_b_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d0;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_psi_ddot;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_psi_ddot = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_psi_ddot),
    &c2_thisId);
  sf_mex_destroy(&c2_psi_ddot);
  *(real_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i56;
  real_T c2_b_inData[4];
  int32_T c2_i57;
  real_T c2_u[4];
  const mxArray *c2_y = NULL;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i56 = 0; c2_i56 < 4; c2_i56++) {
    c2_b_inData[c2_i56] = (*(real_T (*)[4])c2_inData)[c2_i56];
  }

  for (c2_i57 = 0; c2_i57 < 4; c2_i57++) {
    c2_u[c2_i57] = c2_b_inData[c2_i57];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i58;
  real_T c2_b_inData[14];
  int32_T c2_i59;
  real_T c2_u[14];
  const mxArray *c2_y = NULL;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i58 = 0; c2_i58 < 14; c2_i58++) {
    c2_b_inData[c2_i58] = (*(real_T (*)[14])c2_inData)[c2_i58];
  }

  for (c2_i59 = 0; c2_i59 < 14; c2_i59++) {
    c2_u[c2_i59] = c2_b_inData[c2_i59];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 14), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i60;
  real_T c2_b_inData[12];
  int32_T c2_i61;
  real_T c2_u[12];
  const mxArray *c2_y = NULL;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i60 = 0; c2_i60 < 12; c2_i60++) {
    c2_b_inData[c2_i60] = (*(real_T (*)[12])c2_inData)[c2_i60];
  }

  for (c2_i61 = 0; c2_i61 < 12; c2_i61++) {
    c2_u[c2_i61] = c2_b_inData[c2_i61];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 12), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i62;
  real_T c2_b_inData[3];
  int32_T c2_i63;
  real_T c2_u[3];
  const mxArray *c2_y = NULL;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i62 = 0; c2_i62 < 3; c2_i62++) {
    c2_b_inData[c2_i62] = (*(real_T (*)[3])c2_inData)[c2_i62];
  }

  for (c2_i63 = 0; c2_i63 < 3; c2_i63++) {
    c2_u[c2_i63] = c2_b_inData[c2_i63];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static void c2_c_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[3])
{
  real_T c2_dv1[3];
  int32_T c2_i64;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c2_i64 = 0; c2_i64 < 3; c2_i64++) {
    c2_y[c2_i64] = c2_dv1[c2_i64];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_ddot;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[3];
  int32_T c2_i65;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_ddot = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_ddot), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_ddot);
  for (c2_i65 = 0; c2_i65 < 3; c2_i65++) {
    (*(real_T (*)[3])c2_outData)[c2_i65] = c2_y[c2_i65];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i66;
  int32_T c2_i67;
  int32_T c2_i68;
  real_T c2_b_inData[9];
  int32_T c2_i69;
  int32_T c2_i70;
  int32_T c2_i71;
  real_T c2_u[9];
  const mxArray *c2_y = NULL;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_i66 = 0;
  for (c2_i67 = 0; c2_i67 < 3; c2_i67++) {
    for (c2_i68 = 0; c2_i68 < 3; c2_i68++) {
      c2_b_inData[c2_i68 + c2_i66] = (*(real_T (*)[9])c2_inData)[c2_i68 + c2_i66];
    }

    c2_i66 += 3;
  }

  c2_i69 = 0;
  for (c2_i70 = 0; c2_i70 < 3; c2_i70++) {
    for (c2_i71 = 0; c2_i71 < 3; c2_i71++) {
      c2_u[c2_i71 + c2_i69] = c2_b_inData[c2_i71 + c2_i69];
    }

    c2_i69 += 3;
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static void c2_d_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_y[9])
{
  real_T c2_dv2[9];
  int32_T c2_i72;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv2, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c2_i72 = 0; c2_i72 < 9; c2_i72++) {
    c2_y[c2_i72] = c2_dv2[c2_i72];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_T_e2E;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[9];
  int32_T c2_i73;
  int32_T c2_i74;
  int32_T c2_i75;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_T_e2E = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_T_e2E), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_T_e2E);
  c2_i73 = 0;
  for (c2_i74 = 0; c2_i74 < 3; c2_i74++) {
    for (c2_i75 = 0; c2_i75 < 3; c2_i75++) {
      (*(real_T (*)[9])c2_outData)[c2_i75 + c2_i73] = c2_y[c2_i75 + c2_i73];
    }

    c2_i73 += 3;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

const mxArray *sf_c2_quad_dynamics_sim_e_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo;
  c2_ResolvedFunctionInfo c2_info[63];
  const mxArray *c2_m0 = NULL;
  int32_T c2_i76;
  c2_ResolvedFunctionInfo *c2_r0;
  c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  c2_info_helper(c2_info);
  sf_mex_assign(&c2_m0, sf_mex_createstruct("nameCaptureInfo", 1, 63), FALSE);
  for (c2_i76 = 0; c2_i76 < 63; c2_i76++) {
    c2_r0 = &c2_info[c2_i76];
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->context)), "context", "nameCaptureInfo",
                    c2_i76);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c2_r0->name)), "name", "nameCaptureInfo", c2_i76);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c2_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c2_i76);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->resolved)), "resolved", "nameCaptureInfo",
                    c2_i76);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c2_i76);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c2_i76);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c2_i76);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c2_i76);
  }

  sf_mex_assign(&c2_nameCaptureInfo, c2_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c2_nameCaptureInfo);
  return c2_nameCaptureInfo;
}

static void c2_info_helper(c2_ResolvedFunctionInfo c2_info[63])
{
  c2_info[0].context = "";
  c2_info[0].name = "sin";
  c2_info[0].dominantType = "double";
  c2_info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c2_info[0].fileTimeLo = 1343862786U;
  c2_info[0].fileTimeHi = 0U;
  c2_info[0].mFileTimeLo = 0U;
  c2_info[0].mFileTimeHi = 0U;
  c2_info[1].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c2_info[1].name = "eml_scalar_sin";
  c2_info[1].dominantType = "double";
  c2_info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  c2_info[1].fileTimeLo = 1286851136U;
  c2_info[1].fileTimeHi = 0U;
  c2_info[1].mFileTimeLo = 0U;
  c2_info[1].mFileTimeHi = 0U;
  c2_info[2].context = "";
  c2_info[2].name = "mpower";
  c2_info[2].dominantType = "double";
  c2_info[2].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c2_info[2].fileTimeLo = 1286851242U;
  c2_info[2].fileTimeHi = 0U;
  c2_info[2].mFileTimeLo = 0U;
  c2_info[2].mFileTimeHi = 0U;
  c2_info[3].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c2_info[3].name = "power";
  c2_info[3].dominantType = "double";
  c2_info[3].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  c2_info[3].fileTimeLo = 1348224330U;
  c2_info[3].fileTimeHi = 0U;
  c2_info[3].mFileTimeLo = 0U;
  c2_info[3].mFileTimeHi = 0U;
  c2_info[4].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  c2_info[4].name = "eml_scalar_eg";
  c2_info[4].dominantType = "double";
  c2_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[4].fileTimeLo = 1286851196U;
  c2_info[4].fileTimeHi = 0U;
  c2_info[4].mFileTimeLo = 0U;
  c2_info[4].mFileTimeHi = 0U;
  c2_info[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  c2_info[5].name = "eml_scalexp_alloc";
  c2_info[5].dominantType = "double";
  c2_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c2_info[5].fileTimeLo = 1352457260U;
  c2_info[5].fileTimeHi = 0U;
  c2_info[5].mFileTimeLo = 0U;
  c2_info[5].mFileTimeHi = 0U;
  c2_info[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  c2_info[6].name = "floor";
  c2_info[6].dominantType = "double";
  c2_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c2_info[6].fileTimeLo = 1343862780U;
  c2_info[6].fileTimeHi = 0U;
  c2_info[6].mFileTimeLo = 0U;
  c2_info[6].mFileTimeHi = 0U;
  c2_info[7].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c2_info[7].name = "eml_scalar_floor";
  c2_info[7].dominantType = "double";
  c2_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  c2_info[7].fileTimeLo = 1286851126U;
  c2_info[7].fileTimeHi = 0U;
  c2_info[7].mFileTimeLo = 0U;
  c2_info[7].mFileTimeHi = 0U;
  c2_info[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  c2_info[8].name = "eml_scalar_eg";
  c2_info[8].dominantType = "double";
  c2_info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[8].fileTimeLo = 1286851196U;
  c2_info[8].fileTimeHi = 0U;
  c2_info[8].mFileTimeLo = 0U;
  c2_info[8].mFileTimeHi = 0U;
  c2_info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  c2_info[9].name = "mtimes";
  c2_info[9].dominantType = "double";
  c2_info[9].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[9].fileTimeLo = 1289552092U;
  c2_info[9].fileTimeHi = 0U;
  c2_info[9].mFileTimeLo = 0U;
  c2_info[9].mFileTimeHi = 0U;
  c2_info[10].context = "";
  c2_info[10].name = "mtimes";
  c2_info[10].dominantType = "double";
  c2_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[10].fileTimeLo = 1289552092U;
  c2_info[10].fileTimeHi = 0U;
  c2_info[10].mFileTimeLo = 0U;
  c2_info[10].mFileTimeHi = 0U;
  c2_info[11].context = "";
  c2_info[11].name = "cos";
  c2_info[11].dominantType = "double";
  c2_info[11].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c2_info[11].fileTimeLo = 1343862772U;
  c2_info[11].fileTimeHi = 0U;
  c2_info[11].mFileTimeLo = 0U;
  c2_info[11].mFileTimeHi = 0U;
  c2_info[12].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c2_info[12].name = "eml_scalar_cos";
  c2_info[12].dominantType = "double";
  c2_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  c2_info[12].fileTimeLo = 1286851122U;
  c2_info[12].fileTimeHi = 0U;
  c2_info[12].mFileTimeLo = 0U;
  c2_info[12].mFileTimeHi = 0U;
  c2_info[13].context = "";
  c2_info[13].name = "sqrt";
  c2_info[13].dominantType = "double";
  c2_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[13].fileTimeLo = 1343862786U;
  c2_info[13].fileTimeHi = 0U;
  c2_info[13].mFileTimeLo = 0U;
  c2_info[13].mFileTimeHi = 0U;
  c2_info[14].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[14].name = "eml_error";
  c2_info[14].dominantType = "char";
  c2_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  c2_info[14].fileTimeLo = 1343862758U;
  c2_info[14].fileTimeHi = 0U;
  c2_info[14].mFileTimeLo = 0U;
  c2_info[14].mFileTimeHi = 0U;
  c2_info[15].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[15].name = "eml_scalar_sqrt";
  c2_info[15].dominantType = "double";
  c2_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  c2_info[15].fileTimeLo = 1286851138U;
  c2_info[15].fileTimeHi = 0U;
  c2_info[15].mFileTimeLo = 0U;
  c2_info[15].mFileTimeHi = 0U;
  c2_info[16].context = "";
  c2_info[16].name = "mrdivide";
  c2_info[16].dominantType = "double";
  c2_info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c2_info[16].fileTimeLo = 1357983948U;
  c2_info[16].fileTimeHi = 0U;
  c2_info[16].mFileTimeLo = 1319762366U;
  c2_info[16].mFileTimeHi = 0U;
  c2_info[17].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c2_info[17].name = "rdivide";
  c2_info[17].dominantType = "double";
  c2_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[17].fileTimeLo = 1346542788U;
  c2_info[17].fileTimeHi = 0U;
  c2_info[17].mFileTimeLo = 0U;
  c2_info[17].mFileTimeHi = 0U;
  c2_info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[18].name = "eml_scalexp_compatible";
  c2_info[18].dominantType = "double";
  c2_info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  c2_info[18].fileTimeLo = 1286851196U;
  c2_info[18].fileTimeHi = 0U;
  c2_info[18].mFileTimeLo = 0U;
  c2_info[18].mFileTimeHi = 0U;
  c2_info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[19].name = "eml_div";
  c2_info[19].dominantType = "double";
  c2_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c2_info[19].fileTimeLo = 1313380210U;
  c2_info[19].fileTimeHi = 0U;
  c2_info[19].mFileTimeLo = 0U;
  c2_info[19].mFileTimeHi = 0U;
  c2_info[20].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[20].name = "eml_index_class";
  c2_info[20].dominantType = "";
  c2_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[20].fileTimeLo = 1323202978U;
  c2_info[20].fileTimeHi = 0U;
  c2_info[20].mFileTimeLo = 0U;
  c2_info[20].mFileTimeHi = 0U;
  c2_info[21].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[21].name = "eml_scalar_eg";
  c2_info[21].dominantType = "double";
  c2_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[21].fileTimeLo = 1286851196U;
  c2_info[21].fileTimeHi = 0U;
  c2_info[21].mFileTimeLo = 0U;
  c2_info[21].mFileTimeHi = 0U;
  c2_info[22].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[22].name = "eml_xgemm";
  c2_info[22].dominantType = "char";
  c2_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c2_info[22].fileTimeLo = 1299109172U;
  c2_info[22].fileTimeHi = 0U;
  c2_info[22].mFileTimeLo = 0U;
  c2_info[22].mFileTimeHi = 0U;
  c2_info[23].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c2_info[23].name = "eml_blas_inline";
  c2_info[23].dominantType = "";
  c2_info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[23].fileTimeLo = 1299109168U;
  c2_info[23].fileTimeHi = 0U;
  c2_info[23].mFileTimeLo = 0U;
  c2_info[23].mFileTimeHi = 0U;
  c2_info[24].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold";
  c2_info[24].name = "mtimes";
  c2_info[24].dominantType = "double";
  c2_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[24].fileTimeLo = 1289552092U;
  c2_info[24].fileTimeHi = 0U;
  c2_info[24].mFileTimeLo = 0U;
  c2_info[24].mFileTimeHi = 0U;
  c2_info[25].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c2_info[25].name = "eml_index_class";
  c2_info[25].dominantType = "";
  c2_info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[25].fileTimeLo = 1323202978U;
  c2_info[25].fileTimeHi = 0U;
  c2_info[25].mFileTimeLo = 0U;
  c2_info[25].mFileTimeHi = 0U;
  c2_info[26].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c2_info[26].name = "eml_scalar_eg";
  c2_info[26].dominantType = "double";
  c2_info[26].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[26].fileTimeLo = 1286851196U;
  c2_info[26].fileTimeHi = 0U;
  c2_info[26].mFileTimeLo = 0U;
  c2_info[26].mFileTimeHi = 0U;
  c2_info[27].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c2_info[27].name = "eml_refblas_xgemm";
  c2_info[27].dominantType = "char";
  c2_info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c2_info[27].fileTimeLo = 1299109174U;
  c2_info[27].fileTimeHi = 0U;
  c2_info[27].mFileTimeLo = 0U;
  c2_info[27].mFileTimeHi = 0U;
  c2_info[28].context = "";
  c2_info[28].name = "dot";
  c2_info[28].dominantType = "double";
  c2_info[28].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/dot.m";
  c2_info[28].fileTimeLo = 1313380216U;
  c2_info[28].fileTimeHi = 0U;
  c2_info[28].mFileTimeLo = 0U;
  c2_info[28].mFileTimeHi = 0U;
  c2_info[29].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/dot.m";
  c2_info[29].name = "isequal";
  c2_info[29].dominantType = "double";
  c2_info[29].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c2_info[29].fileTimeLo = 1286851158U;
  c2_info[29].fileTimeHi = 0U;
  c2_info[29].mFileTimeLo = 0U;
  c2_info[29].mFileTimeHi = 0U;
  c2_info[30].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c2_info[30].name = "eml_isequal_core";
  c2_info[30].dominantType = "double";
  c2_info[30].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  c2_info[30].fileTimeLo = 1286851186U;
  c2_info[30].fileTimeHi = 0U;
  c2_info[30].mFileTimeLo = 0U;
  c2_info[30].mFileTimeHi = 0U;
  c2_info[31].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar";
  c2_info[31].name = "isnan";
  c2_info[31].dominantType = "double";
  c2_info[31].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  c2_info[31].fileTimeLo = 1286851160U;
  c2_info[31].fileTimeHi = 0U;
  c2_info[31].mFileTimeLo = 0U;
  c2_info[31].mFileTimeHi = 0U;
  c2_info[32].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/dot.m";
  c2_info[32].name = "eml_index_class";
  c2_info[32].dominantType = "";
  c2_info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[32].fileTimeLo = 1323202978U;
  c2_info[32].fileTimeHi = 0U;
  c2_info[32].mFileTimeLo = 0U;
  c2_info[32].mFileTimeHi = 0U;
  c2_info[33].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/dot.m";
  c2_info[33].name = "eml_scalar_eg";
  c2_info[33].dominantType = "double";
  c2_info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[33].fileTimeLo = 1286851196U;
  c2_info[33].fileTimeHi = 0U;
  c2_info[33].mFileTimeLo = 0U;
  c2_info[33].mFileTimeHi = 0U;
  c2_info[34].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/dot.m!vdot";
  c2_info[34].name = "eml_xdotc";
  c2_info[34].dominantType = "double";
  c2_info[34].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m";
  c2_info[34].fileTimeLo = 1299109172U;
  c2_info[34].fileTimeHi = 0U;
  c2_info[34].mFileTimeLo = 0U;
  c2_info[34].mFileTimeHi = 0U;
  c2_info[35].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m";
  c2_info[35].name = "eml_blas_inline";
  c2_info[35].dominantType = "";
  c2_info[35].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[35].fileTimeLo = 1299109168U;
  c2_info[35].fileTimeHi = 0U;
  c2_info[35].mFileTimeLo = 0U;
  c2_info[35].mFileTimeHi = 0U;
  c2_info[36].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m";
  c2_info[36].name = "eml_xdot";
  c2_info[36].dominantType = "double";
  c2_info[36].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m";
  c2_info[36].fileTimeLo = 1299109172U;
  c2_info[36].fileTimeHi = 0U;
  c2_info[36].mFileTimeLo = 0U;
  c2_info[36].mFileTimeHi = 0U;
  c2_info[37].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m";
  c2_info[37].name = "eml_blas_inline";
  c2_info[37].dominantType = "";
  c2_info[37].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[37].fileTimeLo = 1299109168U;
  c2_info[37].fileTimeHi = 0U;
  c2_info[37].mFileTimeLo = 0U;
  c2_info[37].mFileTimeHi = 0U;
  c2_info[38].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m!below_threshold";
  c2_info[38].name = "length";
  c2_info[38].dominantType = "double";
  c2_info[38].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  c2_info[38].fileTimeLo = 1303178606U;
  c2_info[38].fileTimeHi = 0U;
  c2_info[38].mFileTimeLo = 0U;
  c2_info[38].mFileTimeHi = 0U;
  c2_info[39].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m";
  c2_info[39].name = "eml_index_class";
  c2_info[39].dominantType = "";
  c2_info[39].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[39].fileTimeLo = 1323202978U;
  c2_info[39].fileTimeHi = 0U;
  c2_info[39].mFileTimeLo = 0U;
  c2_info[39].mFileTimeHi = 0U;
  c2_info[40].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m";
  c2_info[40].name = "eml_refblas_xdot";
  c2_info[40].dominantType = "double";
  c2_info[40].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c2_info[40].fileTimeLo = 1299109172U;
  c2_info[40].fileTimeHi = 0U;
  c2_info[40].mFileTimeLo = 0U;
  c2_info[40].mFileTimeHi = 0U;
  c2_info[41].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c2_info[41].name = "eml_refblas_xdotx";
  c2_info[41].dominantType = "char";
  c2_info[41].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[41].fileTimeLo = 1299109174U;
  c2_info[41].fileTimeHi = 0U;
  c2_info[41].mFileTimeLo = 0U;
  c2_info[41].mFileTimeHi = 0U;
  c2_info[42].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[42].name = "eml_scalar_eg";
  c2_info[42].dominantType = "double";
  c2_info[42].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[42].fileTimeLo = 1286851196U;
  c2_info[42].fileTimeHi = 0U;
  c2_info[42].mFileTimeLo = 0U;
  c2_info[42].mFileTimeHi = 0U;
  c2_info[43].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[43].name = "eml_index_class";
  c2_info[43].dominantType = "";
  c2_info[43].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[43].fileTimeLo = 1323202978U;
  c2_info[43].fileTimeHi = 0U;
  c2_info[43].mFileTimeLo = 0U;
  c2_info[43].mFileTimeHi = 0U;
  c2_info[44].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[44].name = "eml_int_forloop_overflow_check";
  c2_info[44].dominantType = "";
  c2_info[44].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[44].fileTimeLo = 1346542740U;
  c2_info[44].fileTimeHi = 0U;
  c2_info[44].mFileTimeLo = 0U;
  c2_info[44].mFileTimeHi = 0U;
  c2_info[45].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  c2_info[45].name = "intmax";
  c2_info[45].dominantType = "char";
  c2_info[45].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[45].fileTimeLo = 1311287716U;
  c2_info[45].fileTimeHi = 0U;
  c2_info[45].mFileTimeLo = 0U;
  c2_info[45].mFileTimeHi = 0U;
  c2_info[46].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[46].name = "eml_index_minus";
  c2_info[46].dominantType = "double";
  c2_info[46].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[46].fileTimeLo = 1286851178U;
  c2_info[46].fileTimeHi = 0U;
  c2_info[46].mFileTimeLo = 0U;
  c2_info[46].mFileTimeHi = 0U;
  c2_info[47].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[47].name = "eml_index_class";
  c2_info[47].dominantType = "";
  c2_info[47].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[47].fileTimeLo = 1323202978U;
  c2_info[47].fileTimeHi = 0U;
  c2_info[47].mFileTimeLo = 0U;
  c2_info[47].mFileTimeHi = 0U;
  c2_info[48].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[48].name = "eml_index_times";
  c2_info[48].dominantType = "coder.internal.indexInt";
  c2_info[48].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c2_info[48].fileTimeLo = 1286851180U;
  c2_info[48].fileTimeHi = 0U;
  c2_info[48].mFileTimeLo = 0U;
  c2_info[48].mFileTimeHi = 0U;
  c2_info[49].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c2_info[49].name = "eml_index_class";
  c2_info[49].dominantType = "";
  c2_info[49].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[49].fileTimeLo = 1323202978U;
  c2_info[49].fileTimeHi = 0U;
  c2_info[49].mFileTimeLo = 0U;
  c2_info[49].mFileTimeHi = 0U;
  c2_info[50].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[50].name = "eml_index_plus";
  c2_info[50].dominantType = "coder.internal.indexInt";
  c2_info[50].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[50].fileTimeLo = 1286851178U;
  c2_info[50].fileTimeHi = 0U;
  c2_info[50].mFileTimeLo = 0U;
  c2_info[50].mFileTimeHi = 0U;
  c2_info[51].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[51].name = "eml_index_class";
  c2_info[51].dominantType = "";
  c2_info[51].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[51].fileTimeLo = 1323202978U;
  c2_info[51].fileTimeHi = 0U;
  c2_info[51].mFileTimeLo = 0U;
  c2_info[51].mFileTimeHi = 0U;
  c2_info[52].context = "";
  c2_info[52].name = "mldivide";
  c2_info[52].dominantType = "double";
  c2_info[52].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c2_info[52].fileTimeLo = 1357983948U;
  c2_info[52].fileTimeHi = 0U;
  c2_info[52].mFileTimeLo = 1319762366U;
  c2_info[52].mFileTimeHi = 0U;
  c2_info[53].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c2_info[53].name = "eml_lusolve";
  c2_info[53].dominantType = "double";
  c2_info[53].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c2_info[53].fileTimeLo = 1309483596U;
  c2_info[53].fileTimeHi = 0U;
  c2_info[53].mFileTimeLo = 0U;
  c2_info[53].mFileTimeHi = 0U;
  c2_info[54].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c2_info[54].name = "eml_index_class";
  c2_info[54].dominantType = "";
  c2_info[54].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[54].fileTimeLo = 1323202978U;
  c2_info[54].fileTimeHi = 0U;
  c2_info[54].mFileTimeLo = 0U;
  c2_info[54].mFileTimeHi = 0U;
  c2_info[55].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3";
  c2_info[55].name = "eml_index_class";
  c2_info[55].dominantType = "";
  c2_info[55].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[55].fileTimeLo = 1323202978U;
  c2_info[55].fileTimeHi = 0U;
  c2_info[55].mFileTimeLo = 0U;
  c2_info[55].mFileTimeHi = 0U;
  c2_info[56].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3";
  c2_info[56].name = "eml_xcabs1";
  c2_info[56].dominantType = "double";
  c2_info[56].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c2_info[56].fileTimeLo = 1286851106U;
  c2_info[56].fileTimeHi = 0U;
  c2_info[56].mFileTimeLo = 0U;
  c2_info[56].mFileTimeHi = 0U;
  c2_info[57].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c2_info[57].name = "abs";
  c2_info[57].dominantType = "double";
  c2_info[57].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[57].fileTimeLo = 1343862766U;
  c2_info[57].fileTimeHi = 0U;
  c2_info[57].mFileTimeLo = 0U;
  c2_info[57].mFileTimeHi = 0U;
  c2_info[58].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[58].name = "eml_scalar_abs";
  c2_info[58].dominantType = "double";
  c2_info[58].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c2_info[58].fileTimeLo = 1286851112U;
  c2_info[58].fileTimeHi = 0U;
  c2_info[58].mFileTimeLo = 0U;
  c2_info[58].mFileTimeHi = 0U;
  c2_info[59].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3";
  c2_info[59].name = "eml_div";
  c2_info[59].dominantType = "double";
  c2_info[59].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c2_info[59].fileTimeLo = 1313380210U;
  c2_info[59].fileTimeHi = 0U;
  c2_info[59].mFileTimeLo = 0U;
  c2_info[59].mFileTimeHi = 0U;
  c2_info[60].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3";
  c2_info[60].name = "mtimes";
  c2_info[60].dominantType = "double";
  c2_info[60].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[60].fileTimeLo = 1289552092U;
  c2_info[60].fileTimeHi = 0U;
  c2_info[60].mFileTimeLo = 0U;
  c2_info[60].mFileTimeHi = 0U;
  c2_info[61].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular";
  c2_info[61].name = "eml_warning";
  c2_info[61].dominantType = "char";
  c2_info[61].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m";
  c2_info[61].fileTimeLo = 1286851202U;
  c2_info[61].fileTimeHi = 0U;
  c2_info[61].mFileTimeLo = 0U;
  c2_info[61].mFileTimeHi = 0U;
  c2_info[62].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3";
  c2_info[62].name = "eml_scalar_eg";
  c2_info[62].dominantType = "double";
  c2_info[62].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[62].fileTimeLo = 1286851196U;
  c2_info[62].fileTimeHi = 0U;
  c2_info[62].mFileTimeLo = 0U;
  c2_info[62].mFileTimeHi = 0U;
}

static real_T c2_mpower(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance,
  real_T c2_a)
{
  real_T c2_b_a;
  real_T c2_c_a;
  real_T c2_ak;
  real_T c2_d_a;
  real_T c2_e_a;
  real_T c2_b;
  c2_b_a = c2_a;
  c2_c_a = c2_b_a;
  c2_eml_scalar_eg(chartInstance);
  c2_ak = c2_c_a;
  c2_d_a = c2_ak;
  c2_eml_scalar_eg(chartInstance);
  c2_e_a = c2_d_a;
  c2_b = c2_d_a;
  return c2_e_a * c2_b;
}

static void c2_eml_scalar_eg(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance)
{
}

static void c2_b_eml_scalar_eg(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance)
{
}

static real_T c2_dot(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance,
                     real_T c2_a[3], real_T c2_b[3])
{
  real_T c2_c;
  int32_T c2_k;
  int32_T c2_b_k;
  c2_c_eml_scalar_eg(chartInstance);
  c2_c_eml_scalar_eg(chartInstance);
  c2_c = 0.0;
  for (c2_k = 1; c2_k < 4; c2_k++) {
    c2_b_k = c2_k;
    c2_c += c2_a[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c2_b_k), 1, 3, 1, 0) - 1] * c2_b[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c2_b_k), 1, 3, 1, 0) - 1];
  }

  return c2_c;
}

static void c2_c_eml_scalar_eg(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance)
{
}

static void c2_mldivide(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance,
  real_T c2_A[9], real_T c2_B[3], real_T c2_Y[3])
{
  int32_T c2_i77;
  real_T c2_b_A[9];
  int32_T c2_r1;
  int32_T c2_r2;
  int32_T c2_r3;
  real_T c2_x;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_y;
  real_T c2_d_x;
  real_T c2_e_x;
  real_T c2_b_y;
  real_T c2_maxval;
  real_T c2_f_x;
  real_T c2_g_x;
  real_T c2_h_x;
  real_T c2_c_y;
  real_T c2_i_x;
  real_T c2_j_x;
  real_T c2_d_y;
  real_T c2_a21;
  real_T c2_k_x;
  real_T c2_l_x;
  real_T c2_m_x;
  real_T c2_e_y;
  real_T c2_n_x;
  real_T c2_o_x;
  real_T c2_f_y;
  real_T c2_d;
  real_T c2_p_x;
  real_T c2_g_y;
  real_T c2_z;
  real_T c2_q_x;
  real_T c2_h_y;
  real_T c2_b_z;
  real_T c2_a;
  real_T c2_b;
  real_T c2_i_y;
  real_T c2_b_a;
  real_T c2_b_b;
  real_T c2_j_y;
  real_T c2_c_a;
  real_T c2_c_b;
  real_T c2_k_y;
  real_T c2_d_a;
  real_T c2_d_b;
  real_T c2_l_y;
  real_T c2_r_x;
  real_T c2_s_x;
  real_T c2_t_x;
  real_T c2_m_y;
  real_T c2_u_x;
  real_T c2_v_x;
  real_T c2_n_y;
  real_T c2_b_d;
  real_T c2_w_x;
  real_T c2_x_x;
  real_T c2_y_x;
  real_T c2_o_y;
  real_T c2_ab_x;
  real_T c2_bb_x;
  real_T c2_p_y;
  real_T c2_c_d;
  int32_T c2_rtemp;
  real_T c2_cb_x;
  real_T c2_q_y;
  real_T c2_c_z;
  real_T c2_e_a;
  real_T c2_e_b;
  real_T c2_r_y;
  real_T c2_f_a;
  real_T c2_f_b;
  real_T c2_s_y;
  real_T c2_g_a;
  real_T c2_g_b;
  real_T c2_t_y;
  real_T c2_h_a;
  real_T c2_h_b;
  real_T c2_u_y;
  real_T c2_db_x;
  real_T c2_v_y;
  real_T c2_d_z;
  real_T c2_i_a;
  real_T c2_i_b;
  real_T c2_w_y;
  real_T c2_j_a;
  real_T c2_j_b;
  real_T c2_x_y;
  real_T c2_eb_x;
  real_T c2_y_y;
  real_T c2_e_z;
  real_T c2_k_a;
  real_T c2_k_b;
  real_T c2_ab_y;
  real_T c2_fb_x;
  real_T c2_bb_y;
  real_T c2_f_z;
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  for (c2_i77 = 0; c2_i77 < 9; c2_i77++) {
    c2_b_A[c2_i77] = c2_A[c2_i77];
  }

  c2_r1 = 1;
  c2_r2 = 2;
  c2_r3 = 3;
  c2_x = c2_b_A[0];
  c2_b_x = c2_x;
  c2_c_x = c2_b_x;
  c2_y = muDoubleScalarAbs(c2_c_x);
  c2_d_x = 0.0;
  c2_e_x = c2_d_x;
  c2_b_y = muDoubleScalarAbs(c2_e_x);
  c2_maxval = c2_y + c2_b_y;
  c2_f_x = c2_b_A[1];
  c2_g_x = c2_f_x;
  c2_h_x = c2_g_x;
  c2_c_y = muDoubleScalarAbs(c2_h_x);
  c2_i_x = 0.0;
  c2_j_x = c2_i_x;
  c2_d_y = muDoubleScalarAbs(c2_j_x);
  c2_a21 = c2_c_y + c2_d_y;
  if (c2_a21 > c2_maxval) {
    c2_maxval = c2_a21;
    c2_r1 = 2;
    c2_r2 = 1;
  }

  c2_k_x = c2_b_A[2];
  c2_l_x = c2_k_x;
  c2_m_x = c2_l_x;
  c2_e_y = muDoubleScalarAbs(c2_m_x);
  c2_n_x = 0.0;
  c2_o_x = c2_n_x;
  c2_f_y = muDoubleScalarAbs(c2_o_x);
  c2_d = c2_e_y + c2_f_y;
  if (c2_d > c2_maxval) {
    c2_r1 = 3;
    c2_r2 = 2;
    c2_r3 = 1;
  }

  c2_p_x = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) - 1];
  c2_g_y = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) - 1];
  c2_z = c2_p_x / c2_g_y;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r2), 1, 3, 1, 0) - 1] = c2_z;
  c2_q_x = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) - 1];
  c2_h_y = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) - 1];
  c2_b_z = c2_q_x / c2_h_y;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r3), 1, 3, 1, 0) - 1] = c2_b_z;
  c2_a = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) - 1];
  c2_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) + 2];
  c2_i_y = c2_a * c2_b;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r2), 1, 3, 1, 0) + 2] = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c2_r2), 1, 3, 1, 0) + 2] - c2_i_y;
  c2_b_a = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) - 1];
  c2_b_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) + 2];
  c2_j_y = c2_b_a * c2_b_b;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r3), 1, 3, 1, 0) + 2] = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c2_r3), 1, 3, 1, 0) + 2] - c2_j_y;
  c2_c_a = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) - 1];
  c2_c_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) + 5];
  c2_k_y = c2_c_a * c2_c_b;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r2), 1, 3, 1, 0) + 5] = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c2_r2), 1, 3, 1, 0) + 5] - c2_k_y;
  c2_d_a = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) - 1];
  c2_d_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) + 5];
  c2_l_y = c2_d_a * c2_d_b;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r3), 1, 3, 1, 0) + 5] = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c2_r3), 1, 3, 1, 0) + 5] - c2_l_y;
  c2_r_x = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) + 2];
  c2_s_x = c2_r_x;
  c2_t_x = c2_s_x;
  c2_m_y = muDoubleScalarAbs(c2_t_x);
  c2_u_x = 0.0;
  c2_v_x = c2_u_x;
  c2_n_y = muDoubleScalarAbs(c2_v_x);
  c2_b_d = c2_m_y + c2_n_y;
  c2_w_x = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) + 2];
  c2_x_x = c2_w_x;
  c2_y_x = c2_x_x;
  c2_o_y = muDoubleScalarAbs(c2_y_x);
  c2_ab_x = 0.0;
  c2_bb_x = c2_ab_x;
  c2_p_y = muDoubleScalarAbs(c2_bb_x);
  c2_c_d = c2_o_y + c2_p_y;
  if (c2_b_d > c2_c_d) {
    c2_rtemp = c2_r2;
    c2_r2 = c2_r3;
    c2_r3 = c2_rtemp;
  }

  c2_cb_x = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c2_r3), 1, 3, 1, 0) + 2];
  c2_q_y = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) + 2];
  c2_c_z = c2_cb_x / c2_q_y;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r3), 1, 3, 1, 0) + 2] = c2_c_z;
  c2_e_a = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) + 2];
  c2_e_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) + 5];
  c2_r_y = c2_e_a * c2_e_b;
  c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c2_r3), 1, 3, 1, 0) + 5] = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c2_r3), 1, 3, 1, 0) + 5] - c2_r_y;
  guard1 = FALSE;
  guard2 = FALSE;
  if (c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c2_r1), 1, 3, 1, 0) - 1] == 0.0) {
    guard2 = TRUE;
  } else if (c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
               "", (real_T)c2_r2), 1, 3, 1, 0) + 2] == 0.0) {
    guard2 = TRUE;
  } else {
    if (c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_r3), 1, 3, 1, 0) + 5] == 0.0) {
      guard1 = TRUE;
    }
  }

  if (guard2 == TRUE) {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    c2_eml_warning(chartInstance);
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_Y[0] = c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) - 1];
  c2_f_a = c2_Y[0];
  c2_f_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) - 1];
  c2_s_y = c2_f_a * c2_f_b;
  c2_Y[1] = c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) - 1] - c2_s_y;
  c2_g_a = c2_Y[0];
  c2_g_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) - 1];
  c2_t_y = c2_g_a * c2_g_b;
  c2_h_a = c2_Y[1];
  c2_h_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) + 2];
  c2_u_y = c2_h_a * c2_h_b;
  c2_Y[2] = (c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
               (real_T)c2_r3), 1, 3, 1, 0) - 1] - c2_t_y) - c2_u_y;
  c2_db_x = c2_Y[2];
  c2_v_y = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r3), 1, 3, 1, 0) + 5];
  c2_d_z = c2_db_x / c2_v_y;
  c2_Y[2] = c2_d_z;
  c2_i_a = c2_Y[2];
  c2_i_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) + 5];
  c2_w_y = c2_i_a * c2_i_b;
  c2_Y[0] -= c2_w_y;
  c2_j_a = c2_Y[2];
  c2_j_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) + 5];
  c2_x_y = c2_j_a * c2_j_b;
  c2_Y[1] -= c2_x_y;
  c2_eb_x = c2_Y[1];
  c2_y_y = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r2), 1, 3, 1, 0) + 2];
  c2_e_z = c2_eb_x / c2_y_y;
  c2_Y[1] = c2_e_z;
  c2_k_a = c2_Y[1];
  c2_k_b = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c2_r1), 1, 3, 1, 0) + 2];
  c2_ab_y = c2_k_a * c2_k_b;
  c2_Y[0] -= c2_ab_y;
  c2_fb_x = c2_Y[0];
  c2_bb_y = c2_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c2_r1), 1, 3, 1, 0) - 1];
  c2_f_z = c2_fb_x / c2_bb_y;
  c2_Y[0] = c2_f_z;
}

static void c2_eml_warning(SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance)
{
  int32_T c2_i78;
  static char_T c2_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c2_u[27];
  const mxArray *c2_y = NULL;
  for (c2_i78 = 0; c2_i78 < 27; c2_i78++) {
    c2_u[c2_i78] = c2_varargin_1[c2_i78];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
    14, c2_y));
}

static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(int32_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static int32_T c2_e_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_y;
  int32_T c2_i79;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_i79, 1, 6, 0U, 0, 0U, 0);
  c2_y = c2_i79;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_sfEvent;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_y;
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)chartInstanceVoid;
  c2_b_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_sfEvent),
    &c2_thisId);
  sf_mex_destroy(&c2_b_sfEvent);
  *(int32_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_f_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_quad_dynamics_sim_e, const
  char_T *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_b_is_active_c2_quad_dynamics_sim_e), &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_quad_dynamics_sim_e);
  return c2_y;
}

static uint8_T c2_g_emlrt_marshallIn(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void init_dsm_address_info(SFc2_quad_dynamics_sim_eInstanceStruct
  *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c2_quad_dynamics_sim_e_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(763428315U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(136547201U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(858423591U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3904539814U);
}

mxArray *sf_c2_quad_dynamics_sim_e_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("8B1WLDJ97Ayo2Weaj5p4AH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(14);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c2_quad_dynamics_sim_e_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c2_quad_dynamics_sim_e(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x7'type','srcId','name','auxInfo'{{M[1],M[8],T\"phi_ddot\",},{M[1],M[10],T\"psi_ddot\",},{M[1],M[9],T\"theta_ddot\",},{M[1],M[5],T\"x_ddot\",},{M[1],M[6],T\"y_ddot\",},{M[1],M[7],T\"z_ddot\",},{M[8],M[0],T\"is_active_c2_quad_dynamics_sim_e\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 7, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_quad_dynamics_sim_e_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
    chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *) ((ChartInfoStruct
      *)(ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _quad_dynamics_sim_eMachineNumber_,
           2,
           1,
           1,
           9,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_quad_dynamics_sim_eMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_quad_dynamics_sim_eMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _quad_dynamics_sim_eMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"states");
          _SFD_SET_DATA_PROPS(1,2,0,1,"x_ddot");
          _SFD_SET_DATA_PROPS(2,2,0,1,"y_ddot");
          _SFD_SET_DATA_PROPS(3,2,0,1,"z_ddot");
          _SFD_SET_DATA_PROPS(4,2,0,1,"phi_ddot");
          _SFD_SET_DATA_PROPS(5,2,0,1,"theta_ddot");
          _SFD_SET_DATA_PROPS(6,2,0,1,"psi_ddot");
          _SFD_SET_DATA_PROPS(7,1,1,0,"params");
          _SFD_SET_DATA_PROPS(8,1,1,0,"w");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,6271);
        _SFD_TRANS_COV_WTS(0,0,0,1,0);
        if (chartAlreadyPresent==0) {
          _SFD_TRANS_COV_MAPS(0,
                              0,NULL,NULL,
                              0,NULL,NULL,
                              1,NULL,NULL,
                              0,NULL,NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 14;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c2_x_ddot;
          real_T *c2_y_ddot;
          real_T *c2_z_ddot;
          real_T *c2_phi_ddot;
          real_T *c2_theta_ddot;
          real_T *c2_psi_ddot;
          real_T (*c2_states)[12];
          real_T (*c2_params)[14];
          real_T (*c2_w)[4];
          c2_w = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 2);
          c2_params = (real_T (*)[14])ssGetInputPortSignal(chartInstance->S, 1);
          c2_psi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
          c2_theta_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
          c2_phi_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
          c2_z_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
          c2_y_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c2_x_ddot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          c2_states = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c2_states);
          _SFD_SET_DATA_VALUE_PTR(1U, c2_x_ddot);
          _SFD_SET_DATA_VALUE_PTR(2U, c2_y_ddot);
          _SFD_SET_DATA_VALUE_PTR(3U, c2_z_ddot);
          _SFD_SET_DATA_VALUE_PTR(4U, c2_phi_ddot);
          _SFD_SET_DATA_VALUE_PTR(5U, c2_theta_ddot);
          _SFD_SET_DATA_VALUE_PTR(6U, c2_psi_ddot);
          _SFD_SET_DATA_VALUE_PTR(7U, *c2_params);
          _SFD_SET_DATA_VALUE_PTR(8U, *c2_w);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _quad_dynamics_sim_eMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "dNG29rqBSV13rh8cyxI74C";
}

static void sf_opaque_initialize_c2_quad_dynamics_sim_e(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc2_quad_dynamics_sim_eInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c2_quad_dynamics_sim_e
    ((SFc2_quad_dynamics_sim_eInstanceStruct*) chartInstanceVar);
  initialize_c2_quad_dynamics_sim_e((SFc2_quad_dynamics_sim_eInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c2_quad_dynamics_sim_e(void *chartInstanceVar)
{
  enable_c2_quad_dynamics_sim_e((SFc2_quad_dynamics_sim_eInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c2_quad_dynamics_sim_e(void *chartInstanceVar)
{
  disable_c2_quad_dynamics_sim_e((SFc2_quad_dynamics_sim_eInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c2_quad_dynamics_sim_e(void *chartInstanceVar)
{
  sf_c2_quad_dynamics_sim_e((SFc2_quad_dynamics_sim_eInstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c2_quad_dynamics_sim_e(SimStruct*
  S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c2_quad_dynamics_sim_e
    ((SFc2_quad_dynamics_sim_eInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_quad_dynamics_sim_e();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c2_quad_dynamics_sim_e(SimStruct* S, const
  mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_quad_dynamics_sim_e();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c2_quad_dynamics_sim_e((SFc2_quad_dynamics_sim_eInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c2_quad_dynamics_sim_e(SimStruct*
  S)
{
  return sf_internal_get_sim_state_c2_quad_dynamics_sim_e(S);
}

static void sf_opaque_set_sim_state_c2_quad_dynamics_sim_e(SimStruct* S, const
  mxArray *st)
{
  sf_internal_set_sim_state_c2_quad_dynamics_sim_e(S, st);
}

static void sf_opaque_terminate_c2_quad_dynamics_sim_e(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_quad_dynamics_sim_eInstanceStruct*) chartInstanceVar
      )->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_quad_dynamics_sim_e_optimization_info();
    }

    finalize_c2_quad_dynamics_sim_e((SFc2_quad_dynamics_sim_eInstanceStruct*)
      chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_quad_dynamics_sim_e((SFc2_quad_dynamics_sim_eInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_quad_dynamics_sim_e(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c2_quad_dynamics_sim_e
      ((SFc2_quad_dynamics_sim_eInstanceStruct*)(((ChartInfoStruct *)
         ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c2_quad_dynamics_sim_e(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_quad_dynamics_sim_e_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,2,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,6);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=6; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2312270953U));
  ssSetChecksum1(S,(387573756U));
  ssSetChecksum2(S,(2170947135U));
  ssSetChecksum3(S,(2587834892U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_quad_dynamics_sim_e(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_quad_dynamics_sim_e(SimStruct *S)
{
  SFc2_quad_dynamics_sim_eInstanceStruct *chartInstance;
  chartInstance = (SFc2_quad_dynamics_sim_eInstanceStruct *)utMalloc(sizeof
    (SFc2_quad_dynamics_sim_eInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc2_quad_dynamics_sim_eInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c2_quad_dynamics_sim_e;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c2_quad_dynamics_sim_e_method_dispatcher(SimStruct *S, int_T method, void
  *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_quad_dynamics_sim_e(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_quad_dynamics_sim_e(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_quad_dynamics_sim_e(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_quad_dynamics_sim_e_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
