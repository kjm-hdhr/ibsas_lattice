#include "sign_agg.h"
#include <stdio.h>

void poly_initialize(poly *p){
	int i;
	for(i=0;i<N;i++){
		p->coeffs[i]=0;
	}
}
void polyvecl_initialize(polyvecl *p){
	int i;
	for(i=0;i<L;i++){
		poly_initialize(&(p->vec[i]));
	}
}
void polyveck_initialize(polyveck *p){
	int i;
	for(i=0;i<K;i++){
		poly_initialize(&(p->vec[i]));
	}
}
void copy_polyvecl(polyvecl* dest, const polyvecl* src){
  int i,j;
  for(i=0;i<L;i++){
    for(j=0;j<N;j++){
      dest->vec[i].coeffs[j]=src->vec[i].coeffs[j];
    }
  }
}

void copy_polyveck(polyveck* dest, const polyveck* src){
  int i,j;
  for(i=0;i<K;i++){
    for(j=0;j<N;j++){
      dest->vec[i].coeffs[j]=src->vec[i].coeffs[j];
    }
  }
}
void poly_multiply(poly *c, int32_t a, poly *b, int p){
	int i;
	for(i=0;i<N;i++){
		c->coeffs[i]=multiply(a,b->coeffs[i],p);
	}
}
void poly_pointwise(poly *c, poly *a, poly *b, int p){
	int i;
	for(i=0;i<N;i++){
		c->coeffs[i]=multiply(a->coeffs[i],b->coeffs[i],p);
	}
}
void polyvecl_pointwise_poly(polyvecl *c, poly *a, polyvecl *b, int p){
	int i;
	for(i=0;i<L;i++){
		poly_pointwise(&(c->vec[i]),a,&(b->vec[i]),p);
	}
}
void polyveck_pointwise_poly(polyveck *c, poly *a, polyveck *b, int p){
	int i;
	for(i=0;i<K;i++){
		poly_pointwise(&(c->vec[i]),a,&(b->vec[i]),p);
	}
}
void polyvec_pointwise(polyveck *pk, polyvecl mat[K], polyvecl *pl){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<K;j++){
			pk->vec[j].coeffs[i]=(
				multiply(mat[j].vec[0].coeffs[i],pl->vec[0].coeffs[i],Q)+
				multiply(mat[j].vec[1].coeffs[i],pl->vec[1].coeffs[i],Q)+
				multiply(mat[j].vec[2].coeffs[i],pl->vec[2].coeffs[i],Q)+
				multiply(mat[j].vec[3].coeffs[i],pl->vec[3].coeffs[i],Q)
			)%Q;
		}
	}
}
void poly_add_mod(poly *c, poly *a, poly *b, int p){
	int i;
	for(i=0;i<N;i++){
		c->coeffs[i]=(a->coeffs[i]+b->coeffs[i])%p;
	}
}
void polyvecl_add_mod(polyvecl* c, polyvecl* a, polyvecl* b, int p){
	int i,j;
	for(i=0;i<L;i++){
		poly_add_mod(&(c->vec[i]),&(a->vec[i]),&(b->vec[i]),p);
	}
}
void polyveck_add_mod(polyveck* c, polyveck* a, polyveck* b, int p){
	int i,j;
	for(i=0;i<K;i++){
		poly_add_mod(&(c->vec[i]),&(a->vec[i]),&(b->vec[i]),p);
	}
}
void poly_sub_mod(poly *c, poly *a, poly *b, int p){
	int i;
	for(i=0;i<N;i++){
		c->coeffs[i]=(a->coeffs[i]-b->coeffs[i])%p;
	}
}
void polyvecl_sub_mod(polyvecl* c, polyvecl* a, polyvecl* b, int p){
	int i,j;
	for(i=0;i<L;i++){
		poly_sub_mod(&(c->vec[i]),&(a->vec[i]),&(b->vec[i]),p);
	}
}
void polyveck_sub_mod(polyveck* c, polyveck* a, polyveck* b, int p){
	int i,j;
	for(i=0;i<K;i++){
		poly_sub_mod(&(c->vec[i]),&(a->vec[i]),&(b->vec[i]),p);
	}
}
int compare_polyveck(polyveck *a, polyveck *b){
  for(int i=0;i<K;i++){
    for(int j=0;j<N;j++){
      if(a->vec[i].coeffs[j]!=b->vec[i].coeffs[j]){
        return 0;
	  }
    }
  }
	return 1;
}

int compare_polyvecl(polyvecl *a, polyvecl *b){
  for(int i=0;i<L;i++){
    for(int j=0;j<N;j++){
      if(a->vec[i].coeffs[j]!=b->vec[i].coeffs[j]){
		return 0;
      }
    }
  }
	return 1;
}
int32_t hash_id(uint8_t *id){
	uint8_t idh[4];
	shake256(idh,4,id,SINGLE_ID_LENGTH);
	int32_t* idhash_value=(int32_t*)idh;
	return caddq(*idhash_value%Q);
}
void generate_s1id(polyvecl *s1id, polyvecl mat[K], polyveck *t){
	int32_t mat44[16];
  int32_t mat44_t[16];
  int32_t b[4];
  int32_t inv;
	int i,j;
	for(i=0;i<N;i++){
    for(j=0;j<16;j++){
      mat44[j]=mat[j/4].vec[j%4].coeffs[i];
    }
    for(j=0;j<K;j++){
      b[j]=caddq(t->vec[j].coeffs[i]);
    }
    inv=mat44_det(mat44,Q);
    inv=caddq(inv);
    inv=invmod(inv,Q);
    make_adjugate_mat44(mat44_t,mat44,Q);
    for(j=0;j<K*K;j++){
      mat44_t[j]=caddq(mat44_t[j]);
    }
    for(j=0;j<4;j++){
      s1id->vec[j].coeffs[i]=(
        multiply(mat44_t[j*4+0],b[0],Q)+
        multiply(mat44_t[j*4+1],b[1],Q)+
        multiply(mat44_t[j*4+2],b[2],Q)+
        multiply(mat44_t[j*4+3],b[3],Q))%Q;
      s1id->vec[j].coeffs[i]=multiply(s1id->vec[j].coeffs[i],inv,Q);
    }
  }
}

int crypto_sign_keypair_skid(
	uint8_t *pkid, 
	uint8_t *skid, 
	uint8_t id[SINGLE_ID_LENGTH], 
	uint8_t* sk){
    
	uint8_t rho[SEEDBYTES],tr[SEEDBYTES],key[SEEDBYTES];
	uint8_t trid[SEEDBYTES];
  uint8_t *rhoprime = rho + SEEDBYTES;
	uint8_t idh[4];
  polyvecl mat[K];
  polyvecl s1, idh1, s1id;
  polyveck s2, t,t0, idh2, s2id, pkt, tid;
	uint8_t idhash[CRHBYTES];
  unpack_sk(rho,tr,key,&t0,&s1,&s2,sk);
  polyvec_matrix_expand(mat, rho);

	int32_t idhash_value=hash_id(id);
	shake256(idhash,CRHBYTES,id,SINGLE_ID_LENGTH);
	
	polyvec_pointwise(&t,mat,&s1);
  polyveck_add_mod(&t, &t, &s2,Q);
  polyveck_caddq(&t);

	copy_polyveck(&pkt,&t);
	
	int i,j;
  for(i=0;i<K;i++){
    for(j=0;j<N;j++){
      t.vec[i].coeffs[j]=multiply(t.vec[i].coeffs[j],idhash_value,Q);
      //t.vec[i].coeffs[j]=multiply(t.vec[i].coeffs[j],idhash_value,GAMMA1-1);
    }
  }
	pkid_pack(pkid,rho,&t);
	shake256(trid, SEEDBYTES, pkid, PUBLIC_KEY_AGG_BYTES);

	copy_polyveck(&tid,&t);
	pkid_pack(pkid,rho,&pkt);

  polyveck_uniform_eta(&idh2, idhash, L);
  for(i=0;i<K;i++){
    for(j=0;j<N;j++){
      s2id.vec[i].coeffs[j]=multiply(s2.vec[i].coeffs[j],idh2.vec[i].coeffs[j],Q);
    }
  }
	
  polyveck_sub_mod(&t,&t,&s2id,Q);
	generate_s1id(&s1id, mat, &t);
	
	pkid_pack(pkid,rho,&pkt);
	skid_pack(skid, rho, trid, key, &tid, &s1id, &s2id);
}
int crypto_sign_signatre_skid(
	uint8_t *sig,
	size_t *siglen,
	const uint8_t *m,
	size_t mlen,
	const uint8_t *sk){
	
	unsigned int n;
  uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
	uint8_t c[SEEDBYTES];
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z,s1id;
  polyveck t, t0, s2, w1, w0, h, w,s2id;
  poly cp;
  keccak_state state;

  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + SEEDBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;
	skid_unpack(sk,rho,tr,key,&t,&s1id,&s2id);
  
  /* Compute CRH(tr, msg) */
  shake256_init(&state);
  shake256_absorb(&state, tr, SEEDBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);
  
#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
  shake256(rhoprime, CRHBYTES, key, SEEDBYTES + CRHBYTES);
#endif

  /* Expand matrix and transform vectors */
  //vector A=mat
  polyvec_matrix_expand(mat, rho);
  //polyvecl_ntt(&s1);
  //polyveck_ntt(&s2);
  //polyveck_ntt(&t0);
	int cnt=0;
rej:
  /* Sample intermediate vector y */
  polyvecl_uniform_gamma1(&y, rhoprime, nonce++);
	//polyvecl_uniform_eta(&y, rhoprime, nonce++);
  /* Matrix-vector multiplication */
  z = y;
  //polyvecl_ntt(&z);
	polyvec_pointwise(&w,mat,&y);//Ay
  //polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  //polyveck_reduce(&w1);
  //polyveck_invntt_tomont(&w1);
  
  /* Decompose w and call the random oracle */
  //polyveck_caddq(&w1);
	polyveck_caddq(&w);
  //polyveck_decompose(&w1, &w0, &w1);
  //polyveck_pack_w1(sig, &w1);
	
  polyveck_decompose(&w1, &w0, &w);
  polyveck_pack_w1(sig, &w1);

  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(sig, SEEDBYTES, &state);
	c_pack(c,sig);
	
  //poly_challenge(&cp, sig);
	poly_challenge(&cp, c);
  //poly_ntt(&cp);

  /* Compute z, reject if it reveals secret */
	polyvecl_pointwise_poly(&z, &cp, &s1id,Q);
  //polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
  //polyvecl_invntt_tomont(&z);
  polyvecl_add_mod(&z, &z, &y,Q);
	
  //polyvecl_reduce(&z);
	
  //if(polyvecl_chknorm(&z, GAMMA1 - BETA)){
	//	printf("polyvecl_chknorm(&z, GAMMA1 - BETA) %d\n",cnt++);
  //  goto rej;
	//}

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  //polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
  //polyveck_invntt_tomont(&h);
  //polyveck_sub(&w0, &w0, &h);
  //polyveck_reduce(&w0);
	polyveck cs2id;
	polyveck_pointwise_poly(&cs2id,&cp,&s2id,Q);
	polyveck_sub(&w0,&w0,&cs2id);
	
	//polyveck_pointwise_poly(&h, &cp, &s2id,Q);
	//polyveck_sub(&w0, &w0, &h);
  
	if(polyveck_chknorm(&w0, GAMMA2 - BETA)){
		//printf("polyveck_chknorm(&w0, GAMMA2 - BETA)\n");
    goto rej;
	}
	
  /* Compute hints for w1 */
	/*
  polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
  if(polyveck_chknorm(&h, GAMMA2)){
		printf("polyveck_chknorm(&h, GAMMA2)\n");
    goto rej;
	}

  polyveck_add(&w0, &w0, &h);
  n = polyveck_make_hint(&h, &w0, &w1);
  if(n > OMEGA){
		printf("n > OMEGA\n");
    goto rej;
	}
	*/

  /* Write signature */
	//signature z, c, w
  //pack_sig(sig, sig, &z, &h);
  //*siglen = CRYPTO_BYTES;
	
	sig_pack(sig,c,&z,&w);
	*siglen=SEEDBYTES+POLY_MODQ_LENGTH*L+POLY_MODQ_LENGTH*K;
  return 0;
}

int crypto_sign_signature_agg(
	uint8_t *agg_sig,
	uint8_t *prev_agg_sig,
	uint8_t *sig,
	uint8_t prev_agg_sig_wlen){
	
	uint8_t c_p[SEEDBYTES];
	polyvecl z_p;
	polyveck w_p[prev_agg_sig_wlen+1];
	agg_sig_unpack(prev_agg_sig,c_p,&z_p, w_p,prev_agg_sig_wlen);

	uint8_t c[SEEDBYTES];
	polyvecl z;
	polyveck w;
	sig_unpack(sig,c,&z,&w);
	w_p[prev_agg_sig_wlen]=w;
	polyvecl_add_mod(&z, &z, &z_p,Q);
	int i;
	for(i=0;i<SEEDBYTES;i++){
		c[i]=(c[i]+c_p[i])%Q;
	}
	agg_sig_pack(agg_sig,c,&z,w_p,prev_agg_sig_wlen+1);
}

int crypto_sign_verify_agg(uint8_t *pk, uint8_t *agg_sig, uint8_t *msg, int idlen){
	uint8_t c_agg[SEEDBYTES];
	uint8_t cid[SEEDBYTES]={0},c_tmp[SEEDBYTES];
	polyvecl z_agg;
	polyvecl mat[K];
	polyveck w_agg[idlen];
	polyveck w1,w0,wid,wsum;
	polyveck ctid,ct_tmp;
	polyveck Az;
	poly pc;
	polyveck_initialize(&wsum);
	polyveck_initialize(&ctid);
	int wlen;
	agg_sig_unpack(agg_sig,c_agg,&z_agg, w_agg,idlen);
	polyveck t,tid;
	uint8_t rho[SEEDBYTES];
	pkid_unpack(pk,rho,&t);//rho=A, t=t

	polyvec_matrix_expand(mat, rho);
	polyvec_pointwise(&Az,mat,&z_agg);
	polyveck_caddq(&Az);
	
	int i,j;
	int32_t idhash_value;
	uint8_t pkid[PUBLIC_KEY_AGG_BYTES];
	uint8_t trid[SEEDBYTES];
	keccak_state state;
	uint8_t mu[CRHBYTES];
	uint8_t w1_pack[K*POLYW1_PACKEDBYTES];

	uint8_t id_tmp[SINGLE_ID_LENGTH];
	for(i=1;i<=idlen;i++){
		id_tmp[0]=msg[(i-1)*SINGLE_ID_LENGTH+0];
		id_tmp[1]=msg[(i-1)*SINGLE_ID_LENGTH+1];
		id_tmp[2]=msg[(i-1)*SINGLE_ID_LENGTH+2];
		id_tmp[3]=msg[(i-1)*SINGLE_ID_LENGTH+3];
		idhash_value=hash_id(id_tmp);
		for(j=0;j<K;j++){
			poly_multiply(&(tid.vec[j]),idhash_value,&(t.vec[j]),Q);
		}
		
		pkid_pack(pkid,rho,&tid);
		// tr = CRH(rho || t)
		shake256(trid, SEEDBYTES, pkid, PUBLIC_KEY_AGG_BYTES);

		shake256_init(&state);
  	shake256_absorb(&state, trid, SEEDBYTES);
  	shake256_absorb(&state, msg, i*SINGLE_ID_LENGTH);
  	shake256_finalize(&state);
		
		//mu = CRH(tr||M)
  	shake256_squeeze(mu, CRHBYTES, &state);
		polyveck_decompose(&w1, &w0, &w_agg[i-1]);
		polyveck_pack_w1(w1_pack, &w1);

		shake256_init(&state);
  	shake256_absorb(&state, mu, CRHBYTES);
  	shake256_absorb(&state, w1_pack, K*POLYW1_PACKEDBYTES);
  	shake256_finalize(&state);
		//c = CRH(mu||w1)
  	shake256_squeeze(c_tmp, SEEDBYTES, &state);
		for(j=0;j<SEEDBYTES;j++){
			cid[j]=(cid[j]+c_tmp[j])%Q;
		}
		//ct
		poly_challenge(&pc, c_tmp);
		polyveck_add_mod(&wsum,&wsum,&w_agg[i-1],Q);
		polyveck_pointwise_poly(&ct_tmp,&pc,&tid,Q);
		polyveck_add_mod(&ctid,&ctid,&ct_tmp,Q);
	}
	for(i=0;i<SEEDBYTES;i++){
		if(c_agg[i]!=cid[i]){	
			return 0;
		}
	}
	//Az = tc+w
	polyveck_sub_mod(&wid,&Az,&ctid,Q);
	polyveck_caddq(&wid);
	
	polyveck wid1, wid0;
	polyveck wsum1, wsum0;
	polyveck_decompose(&wid1, &wid0, &wid);
	polyveck_decompose(&wsum1, &wsum0, &wsum);
	return compare_polyveck(&wsum1,&wid1);
	
}

	