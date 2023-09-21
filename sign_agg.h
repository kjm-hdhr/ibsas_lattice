#ifndef SIGN_AGG_H
#define SIGN_AGG_H
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "randombytes.h"
#include "sign.h"
#include "packing.h"
#include "fips202.h"
#include <time.h>
#include "matrix.h"
#include "reduce.h"
#include "packing_agg.h"


#define SINGLE_ID_LENGTH 4

typedef struct secret_key_id{
	uint8_t rho[SEEDBYTES];
	uint8_t tr[SEEDBYTES];
	uint8_t key[SEEDBYTES];
	uint8_t t[POLY_MODQ_LENGTH*K];
	uint8_t s1id[POLY_MODQ_LENGTH*L];
	uint8_t s2id[POLYS2ID_LENGTH*4];
} sec_key_id;

typedef struct sig{
	polyvecl z;
	polyveck w;
} sig;

typedef struct agg_sig{
	polyvecl z;
	polyveck *w;
	int wlen;
} agg_sig;

void poly_initialize(poly *p);
void polyvecl_initialize(polyvecl *p);
void polyveck_initialize(polyveck *p);
void copy_polyvecl(polyvecl* dest, const polyvecl* src);
void copy_polyveck(polyveck* dest, const polyveck* src);
void poly_pointwise(poly *c, poly *a, poly *b, int p);
void polyvecl_pointwise_poly(polyvecl *c, poly *a, polyvecl *b, int p);
void polyveck_pointwise_poly(polyveck *c, poly *a, polyveck *b, int p);
void polyvec_pointwise(polyveck *pk, polyvecl mat[K], polyvecl *pl);
void poly_add_mod(poly *c, poly *a, poly *b, int p);
void polyvecl_add_mod(polyvecl* c, polyvecl* a, polyvecl* b, int p);
void polyveck_add_mod(polyveck* c, polyveck* a, polyveck* b, int p);
void poly_sub_mod(poly *c, poly *a, poly *b, int p);
void polyvecl_sub_mod(polyvecl* c, polyvecl* a, polyvecl* b, int p);
void polyveck_sub_mod(polyveck* c, polyveck* a, polyveck* b, int p);
void poly_multiply(poly *c, int32_t a, poly *b, int p);
int compare_polyvecl(polyvecl *a, polyvecl *b);
int compare_polyveck(polyveck *a, polyveck *b);
void generate_s1id(polyvecl *s1id, polyvecl mat[K], polyveck *t);
int32_t hash_id(uint8_t *id);

void setup_skid();
int crypto_sign_keypair_skid(
	uint8_t *pkid, 
	uint8_t *skid, 
	uint8_t id[SINGLE_ID_LENGTH],
	uint8_t* sk);

int crypto_sign_signatre_skid(
	uint8_t *sig,
	size_t *siglen,
	const uint8_t *m,
	size_t mlen,
	const uint8_t *sk);

int crypto_sign_signature_agg(
	uint8_t *agg_sig,
	uint8_t *prev_agg_sig,
	uint8_t *sig,
	uint8_t prev_agg_sig_mlen);

int crypto_sign_verify_agg(uint8_t *pk, uint8_t *agg_sig, uint8_t *msg, int idlen);
#endif