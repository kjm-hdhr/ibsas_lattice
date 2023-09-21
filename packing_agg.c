#include "packing_agg.h"
#include <stdio.h>
void c_pack(uint8_t *pack, const uint8_t c[SEEDBYTES]){
	int i;
	for(i=0;i<SEEDBYTES;i++){
		pack[i]=c[i];
	}
}
void c_unpack(const uint8_t *pack, uint8_t c[SEEDBYTES]){
	int i;
	for(i=0;i<SEEDBYTES;i++){
		c[i]=pack[i];
	}
}
void poly_pack_modQ(uint8_t *pack, const poly *p){
	int i;
	for(i=0;i<N;i++){
		pack[i*3+0]=p->coeffs[i];
		pack[i*3+1]=p->coeffs[i]>>8;
		pack[i*3+2]=p->coeffs[i]>>16;
		//printf("pack[%x %x %x]=coeffs[%d] %x\n",pack[i*3+0],pack[i*3+1],pack[i*3+2],i,p->coeffs[i]);
	}
}
void poly_unpack_modQ(const uint8_t *pack, poly *p){
	int i;
	for(i=0;i<N;i++){
		p->coeffs[i]=pack[i*3+0];
		p->coeffs[i]|=((int32_t)pack[i*3+1])<<8;
		p->coeffs[i]|=(((int32_t)pack[i*3+2])<<24)>>8;
		//printf("unpack[%x %x %x]=coeffs[%d] %x\n",pack[i*3+0],pack[i*3+1],pack[i*3+2],i,p->coeffs[i]);
	}
}

void polyveck_pack_modQ(uint8_t *pack, const polyveck *k){
	int i;
	for(i=0;i<K;i++){
		poly_pack_modQ(pack,&(k->vec[i]));
		pack+=POLY_MODQ_LENGTH;
	}
}
void polyveck_unpack_modQ(const uint8_t *pack, polyveck *k){
	int i;
	for(i=0;i<K;i++){
		poly_unpack_modQ(pack,&(k->vec[i]));
		pack+=POLY_MODQ_LENGTH;
	}
}
void polyvecl_pack_modQ(uint8_t *pack, const polyvecl *l){
	int i;
	for(i=0;i<L;i++){
		poly_pack_modQ(pack,&(l->vec[i]));
		pack+=POLY_MODQ_LENGTH;
	}
}
void polyvecl_unpack_modQ(const uint8_t *pack, polyvecl *l){
	int i;
	for(i=0;i<L;i++){
		poly_unpack_modQ(pack,&(l->vec[i]));
		pack+=POLY_MODQ_LENGTH;
	}
}

void polyvecz_pack(uint8_t *pack, const polyvecl *z){
	polyvecl_pack_modQ(pack,z);
}
void polyvecz_unpack(const uint8_t *pack, polyvecl *z){
	polyvecl_unpack_modQ(pack,z);
}
void polyvecw_pack(uint8_t *pack, const polyveck *w){
	polyveck_pack_modQ(pack,w);
}
void polyvecw_unpack(const uint8_t *pack, polyveck *w){
	polyveck_unpack_modQ(pack,w);
}
void polyvect_pack(uint8_t *pack, const polyveck *t){
	polyveck_pack_modQ(pack,t);
}
void polyvect_unpack(const uint8_t *pack, polyveck *t){
	polyveck_unpack_modQ(pack,t);
}
void polyvecs1id_pack(uint8_t *pack, const polyvecl *s1id){
	polyvecl_pack_modQ(pack,s1id);
}
void polyvecs1id_unpack(const uint8_t *pack, polyvecl *s1id){
	polyvecl_unpack_modQ(pack,s1id);
}
void polys2id_pack(uint8_t *pack, const poly *s2id){
	uint8_t t[2];
	int i;
	for(i = 0; i < N/2; ++i) {
    t[0] = 4 - s2id->coeffs[2*i+0];
    t[1] = 4 - s2id->coeffs[2*i+1];
    pack[i] = t[0] | (t[1] << 4);
  }
}
void polys2id_unpack(const uint8_t *pack, poly *s2id){
	int i;
	for(i = 0; i < N/2; ++i) {
    s2id->coeffs[2*i+0] = pack[i] & 0x0F;
    s2id->coeffs[2*i+1] = pack[i] >> 4;
    s2id->coeffs[2*i+0] = 4 - s2id->coeffs[2*i+0];
    s2id->coeffs[2*i+1] = 4 - s2id->coeffs[2*i+1];
  }
}
void polyvecs2id_pack(uint8_t *pack, const polyveck *s2id){
	int i;
	for(i=0;i<K;i++){
		polys2id_pack(pack,&(s2id->vec[i]));
		pack+=POLYS2ID_LENGTH;
	}
}
void polyvecs2id_unpack(const uint8_t *pack, polyveck *s2id){
	int i;
	for(i=0;i<K;i++){
		polys2id_unpack(pack,&(s2id->vec[i]));
		pack+=POLYS2ID_LENGTH;
	}
}
void pkid_pack(uint8_t *pkid, uint8_t *rho, polyveck *t){
	int i;
	for(i=0;i<SEEDBYTES;i++){
		pkid[i]=rho[i];
	}
	pkid+=SEEDBYTES;
	polyvect_pack(pkid,t);
}

void pkid_unpack(uint8_t *pkid, uint8_t *rho, polyveck *t){
	int i;
	for(i=0;i<SEEDBYTES;i++){
		rho[i]=pkid[i];
	}
	pkid+=SEEDBYTES;
	polyvect_unpack(pkid,t);
}
void skid_pack(
	uint8_t *skid, 
	uint8_t* rho, uint8_t *tr,uint8_t *key,
	polyveck *t, polyvecl *s1id, polyveck *s2id){
	
	int i;
	for(i=0;i<SEEDBYTES;i++){
		skid[i]=rho[i];
		skid[i+SEEDBYTES]=tr[i];
		skid[i+SEEDBYTES*2]=key[i];
	}
	skid+=SEEDBYTES*3;
	polyvect_pack(skid,t);
	skid+=POLY_MODQ_LENGTH*K;
	polyvecs1id_pack(skid,s1id);
	skid+=POLY_MODQ_LENGTH*L;
	polyvecs2id_pack(skid,s2id);
}
void skid_unpack(const uint8_t *skid, 
	uint8_t* rho, uint8_t *tr,uint8_t *key,
	polyveck *t, polyvecl *s1id, polyveck *s2id){
	
	int i;
	for(i=0;i<SEEDBYTES;i++){
		rho[i]=skid[i];
		tr[i]=skid[i+SEEDBYTES];
		key[i]=skid[i+SEEDBYTES*2];
	}
	skid+=SEEDBYTES*3;
	polyvect_unpack(skid,t);
	skid+=POLY_MODQ_LENGTH*K;
	polyvecs1id_unpack(skid,s1id);
	skid+=POLY_MODQ_LENGTH*L;
	polyvecs2id_unpack(skid,s2id);
}

void sig_pack(uint8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w){
	c_pack(pack,c);
	pack+=SEEDBYTES;
	polyvecz_pack(pack,z);
	pack+=POLY_MODQ_LENGTH*L;
	polyvecw_pack(pack,w);
}
void sig_unpack(uint8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w){
	c_unpack(pack,c);
	pack+=SEEDBYTES;
	polyvecz_unpack(pack,z);
	pack+=POLY_MODQ_LENGTH*L;
	polyvecw_unpack(pack,w);
}
void agg_sig_pack(int8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w, uint8_t wlen){
	c_pack(pack,c);
	pack+=SEEDBYTES;
	polyvecz_pack(pack,z);
	pack+=POLY_MODQ_LENGTH*L;
	int i;
	for(i=0;i<wlen;i++){
		//w+=i;
		polyvecw_pack(pack,&w[i]);
		pack+=POLY_MODQ_LENGTH*K;
	}
}
void agg_sig_unpack(int8_t *pack, uint8_t c[SEEDBYTES], polyvecl *z, polyveck *w, uint8_t wlen){
	c_unpack(pack,c);
	pack+=SEEDBYTES;
	polyvecz_unpack(pack,z);
	pack+=POLY_MODQ_LENGTH*L;
	int i;
	for(i=0;i<wlen;i++){
		//w+=i;
		polyvecw_unpack(pack,&w[i]);
		pack+=POLY_MODQ_LENGTH*K;
	}
}
