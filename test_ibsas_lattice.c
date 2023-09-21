#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "randombytes.h"
#include "sign.h"
#include "packing.h"
#include "fips202.h"
#include <time.h>
#include "matrix.h"
#include "reduce.h"
#include "sign_agg.h"
//#include "pico/stdlib.h"

#define SINGLE_MSGLEN 4
#define NODES 10
#define REPEAT 100

int main(void)
{
  clock_t t;

  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];

  double kg[NODES]={0},sg[NODES]={0},sa[NODES]={0},sv[NODES]={0};
  uint8_t pkid[PUBLIC_KEY_AGG_BYTES];
  uint8_t skid[SECRET_KEY_AGG_BYTES];
  uint8_t current[SINGLE_ID_LENGTH]={10,0,0,0};
  uint8_t sig[EACH_SIG_LENGTH];
  uint8_t *aggsig,*aggsig_prev;
  size_t siglen,msglen;
  int r=0;
  int i,j,v;
  crypto_sign_keypair(pk,sk);

  for(r=0;r<REPEAT;r++){
    printf("repeat %d\n",(r+1));

    // processes at the first node
    memset(pkid,0,PUBLIC_KEY_AGG_BYTES);
    memset(skid,0,SECRET_KEY_AGG_BYTES);
    memset(sig,0,EACH_SIG_LENGTH);
    aggsig=malloc(SEEDBYTES+POLY_MODQ_LENGTH*L+POLY_MODQ_LENGTH*K+1);
    memset(aggsig,0,SEEDBYTES+POLY_MODQ_LENGTH*L+POLY_MODQ_LENGTH*K+1);
    current[3]=1;
    msglen=SINGLE_ID_LENGTH;
    uint8_t *msg=malloc(SINGLE_ID_LENGTH);
    msg[0]=10,msg[1]=0,msg[2]=0,msg[3]=1;
    t=clock();
    crypto_sign_keypair_skid(pkid, skid,current,sk);
    kg[0]+=clock()-t;

    t=clock();
    crypto_sign_signatre_skid(
    sig, &siglen,msg,msglen,skid);
    sg[0]+=clock()-t;

    t=clock();
    crypto_sign_signature_agg(aggsig,aggsig,sig,0);
    sg[0]+=clock()-t;

    t=clock();
    v=crypto_sign_verify_agg(pk,aggsig,msg,1);
    sv[0]+=clock()-t;

    free(msg);
    aggsig_prev=aggsig;

    for(i=1;i<NODES;i++){
      memset(pkid,0,PUBLIC_KEY_AGG_BYTES);
      memset(skid,0,SECRET_KEY_AGG_BYTES);
      memset(sig,0,EACH_SIG_LENGTH);
      aggsig=malloc(SEEDBYTES+POLY_MODQ_LENGTH*L+POLY_MODQ_LENGTH*K*(i+1)+1);
      memset(aggsig,0,SEEDBYTES+POLY_MODQ_LENGTH*L+POLY_MODQ_LENGTH*K*(i+1)+1);
      current[3]=(i+1);
      msglen=SINGLE_ID_LENGTH*(i+1);
      uint8_t *msg=malloc(SINGLE_ID_LENGTH*(i+1));
      for(j=0;j<(i+1);j++){
        msg[(SINGLE_ID_LENGTH*j)+0]=10;
        msg[(SINGLE_ID_LENGTH*j)+1]=0;
        msg[(SINGLE_ID_LENGTH*j)+2]=0;
        msg[(SINGLE_ID_LENGTH*j)+3]=(j+1);
      }
      t=clock();
      crypto_sign_keypair_skid(pkid, skid,current,sk);
      kg[i]+=clock()-t;

      t=clock();
      crypto_sign_signatre_skid(
      sig, &siglen,msg,msglen,skid);
      sg[i]+=clock()-t;

      t=clock();
      crypto_sign_signature_agg(aggsig,aggsig_prev,sig,i);
      sa[i]+=clock()-t;

      t=clock();
      v=crypto_sign_verify_agg(pk,aggsig,msg,(i+1));
      sv[i]+=clock()-t;
      free(msg);
      free(aggsig_prev);
      aggsig_prev=aggsig;
    }
  }
  double kgavg[NODES]={0.0},sgavg[NODES]={0.0},saavg[NODES]={0.0},svavg[NODES]={0.0};
  for(i=0;i<NODES;i++){
    kgavg[i]=(kg[i]/REPEAT)/CLOCKS_PER_SEC;
    sgavg[i]=(sg[i]/REPEAT)/CLOCKS_PER_SEC;
    saavg[i]=(sa[i]/REPEAT)/CLOCKS_PER_SEC;
    svavg[i]=(sv[i]/REPEAT)/CLOCKS_PER_SEC;
    //if((i+1)%10==0)
    printf("i=%d,kg=%f,sg=%f,sa=%f,sg+sa=%f,sv=%f\n",(i+1),kgavg[i],sgavg[i],saavg[i],sgavg[i]+saavg[i],svavg[i]);
  }
  return 0;
}
