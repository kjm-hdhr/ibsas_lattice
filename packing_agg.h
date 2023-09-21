#ifndef PACKING_AGG_H
#define PACKING_AGG_H
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

#define POLY_MODQ_LENGTH 256*3
#define POLYS2ID_LENGTH 128
#define SECRET_KEY_AGG_BYTES SEEDBYTES*3+POLY_MODQ_LENGTH*K+POLY_MODQ_LENGTH*L+POLYS2ID_LENGTH*K
#define PUBLIC_KEY_AGG_BYTES SEEDBYTES+POLY_MODQ_LENGTH*K
#define EACH_SIG_LENGTH SEEDBYTES+POLY_MODQ_LENGTH*L+POLY_MODQ_LENGTH*K

void poly_pack_modQ(uint8_t *pack, const poly *polyt);
void poly_unpack_modQ(const uint8_t *pack, poly *polyt);

void polyveck_pack_modQ(uint8_t *pack, const polyveck *k);
void polyveck_unpack_modQ(const uint8_t *pack, polyveck *k);
void polyvecl_pack_modQ(uint8_t *pack, const polyvecl *l);
void polyvecl_unpack_modQ(const uint8_t *pack, polyvecl *l);

void polyvecz_pack(uint8_t *pack, const polyvecl *z);
void polyvecz_unpack(const uint8_t *pack, polyvecl *z);
void polyvecw_pack(uint8_t *pack, const polyveck *w);
void polyvecw_unpack(const uint8_t *pack, polyveck *w);

void polyvect_pack(uint8_t *pack, const polyveck *t);
void polyvect_unpack(const uint8_t *pack, polyveck *t);
void polyvecs1id_pack(uint8_t *pack, const polyvecl *s1id);
void polyvecs1id_unpack(const uint8_t *pack, polyvecl *s1id);

void polys2id_pack(uint8_t *pack, const poly *s2id);
void polys2id_unpack(const uint8_t *pack, poly *s2id);
void polyvecs2id_pack(uint8_t *pack, const polyveck *s2id);
void polyvecs2id_unpack(const uint8_t *pack, polyveck *s2id);

void sig_pack(uint8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w);
void sig_unpack(uint8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w);
void agg_sig_pack(int8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w, uint8_t wlen);
void agg_sig_unpack(int8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w, uint8_t wlen);
void c_pack(uint8_t *pack, const uint8_t c[SEEDBYTES]);
void c_unpack(const uint8_t *pack, uint8_t c[SEEDBYTES]);

void skid_pack(uint8_t *skid, uint8_t* rho, uint8_t *tr,uint8_t *key,polyveck *t, polyvecl *s1id, polyveck *s2id);
void skid_unpack(const uint8_t *skid, uint8_t* rho, uint8_t *tr,uint8_t *key,polyveck *t, polyvecl *s1id, polyveck *s2id);
void pkid_pack(uint8_t *pkid, uint8_t *rho, polyveck *t);
void pkid_unpack(uint8_t *pkid, uint8_t *rho, polyveck *t);
#endif