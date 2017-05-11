#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>

#define TRUE 1
#define FALSE 0

mpz_t prime;
mpz_t trace;
mpz_t order;
mpz_t r;

mpz_t ECCPara_a;//y^2 = x^3+ax^2+x

struct EFp{
	mpz_t x,y;
	int Inf;
};

struct EFp_XZ{
	mpz_t x,z;
	int Inf;
}

//----------Fp functions--------------------------------------------
void Fp_Pow(mpz_t Ans, mpz_t A, mpz_t B);	//Ans <- A^B
void Fp_Inv(mpz_t Ans, mpz_t A);			//Ans <- A^(-1)
void Fp_Sqrt(mpz_t Ans, mpz_t A);			//ANs <- A^(1/2)
//----------EFp functions--------------------------------------------
void EFp_Init(struct EFp *A);	// vector initialization
void EFp_Copy(struct EFp *A, struct EFp *B);	//A <- B
void EFp_SetInf(struct EFp *A);	//A <- Infinity
void EFp_Clear(struct EFp *A);
void EFp_Rand(struct EFp *A);
void EFp_Show(struct EFp *A);
void EFp_ECA(struct EFp *Ans, struct EFp *P, struct EFp *Q);//ANS=P+Q
void EFp_ECD(struct EFp *Ans, struct EFp *P);//ANS=2*P
void EFp_SCM(struct EFp *Ans, struct EFp *P,mpz_t j);//ANS=[j]P

void EFp_XZ_Init(struct EFp_XZ *A);
void EFp_XZ_Clear(struct EFp_XZ *A);
void EFp_XZ_ECA(struct EFp_XZ *Ans, struct EFp_XZ *P, struct EFp_XZ *Q);//ECA by XZ 
void EFp_XZ_ECD(struct EFp_XZ *Ans, struct EFp_XZ *P);//ECD by XZ

//-------------Fp functions--------------------------------
void Fp_Inv(mpz_t Ans, mpz_t A){
	mpz_t Ans_Copy,j;
	mpz_init(Ans_Copy);
	mpz_init(j);
	mpz_set(Ans_Copy, A);
	mpz_sub_ui(j,prime,2);
	int i,r;
	r=(int)mpz_sizeinbase(j,2);
	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);
			mpz_mod(Ans_Copy,Ans_Copy,prime);
			mpz_mul(Ans_Copy,Ans_Copy,A);
		}else{
		mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);
		mpz_mod(Ans_Copy,Ans_Copy,prime);
		}
	}
	mpz_set(Ans,Ans_Copy);

	mpz_clear(Ans_Copy);
	mpz_clear(j);
}

void Fp_Pow(mpz_t Ans, mpz_t A, mpz_t B){
	if(mpz_cmp_ui(B,0)==0){
		mpz_set_ui(Ans,1);
		return;
	}
	mpz_t Ans_Copy;
	int i;
	int r;//bit

	mpz_init(Ans_Copy);
	r= (int)mpz_sizeinbase(B,2);
	mpz_set(Ans_Copy,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);	
			mpz_mul(Ans_Copy,Ans_Copy,A);			//(A*2)*A
			mpz_mod(Ans_Copy,Ans_Copy,prime);
		}else{
			mpz_mul(Ans_Copy,Ans_Copy,Ans_Copy);	//A*2
			mpz_mod(Ans_Copy,Ans_Copy,prime);
		}
	}
	mpz_mod(Ans_Copy,Ans_Copy,prime);
	mpz_set(Ans,Ans_Copy);

	mpz_clear(Ans_Copy);
}

void Fp_Sqrt(mpz_t Ans, mpz_t A){
	if(mpz_legendre(A,prime)!=1){
		printf("No Answer¥n");
		return;
	}

	mpz_t b,e,n,q,r,t,x,y,tmp1,tmp2,base;
	unsigned int m,Rand;
	gmp_randstate_t state;
	gmp_randinit_default(state);

	mpz_init(b);
	mpz_init(e);
	mpz_init(n);
	mpz_init(q);
	mpz_init(r);
	mpz_init(t);
	mpz_init(x);
	mpz_init(y);
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(base);
	
	mpz_set_ui(base,2);
	mpz_set(n,A);
	while(mpz_legendre(n,prime)!=-1){
		Rand=rand();
		gmp_randseed_ui(state,Rand);
		mpz_urandomm(n, state, prime);
	}
	mpz_set_ui(n,21);
	mpz_set_ui(e,0);
	mpz_sub_ui(q,prime,1);
	
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	mpz_powm(y,n,q,prime);
	mpz_set(r,e);
	mpz_sub_ui(q,q,1);
	mpz_div_ui(q,q,2);
	mpz_powm(x,A,q,prime);
	mpz_set(tmp1,x);
	mpz_mul(x,x,A);
	mpz_mul(b,x,tmp1);
	mpz_mod(b,b,prime);

	while(mpz_cmp_ui(b,1)){
		m=-1;
		mpz_init(tmp2);
		while(mpz_cmp_ui(tmp2,1)){
			m++;
			mpz_ui_pow_ui(tmp1,2,m);
			mpz_powm(tmp2,b,tmp1,prime);
		}
		mpz_sub_ui(tmp1,r,m);
		mpz_sub_ui(tmp1,tmp1,1);
		mpz_powm(tmp2,base,tmp1,prime);
		mpz_powm(t,y,tmp2,prime);
		mpz_powm_ui(y,t,2,prime);
		mpz_set_ui(r,m);
		mpz_mul(x,x,t);
		mpz_mod(x,x,prime);
		mpz_mul(b,b,y);
		mpz_mod(b,b,prime);
	}
	
	mpz_set(Ans,x);
	mpz_mod(Ans,Ans,prime);
	
	mpz_clear(b);
	mpz_clear(e);
	mpz_clear(n);
	mpz_clear(q);
	mpz_clear(r);
	mpz_clear(t);
	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(tmp2);
	mpz_clear(tmp1);
	mpz_clear(base);

}

//-------------EFp functions--------------------------------
void EFp_Init(struct EFp *A){
	mpz_init(A->x);
	mpz_init(A->y);
	A->Inf=FALSE;
}

void EFp_Copy(struct EFp *A, struct EFp *B){
	mpz_set(A->x,B->x);
	mpz_set(A->y,B->y);
	A->Inf=B->Inf;
}

void EFp_SetInf(struct EFp *A){
	mpz_set_ui(A->x,0);
	mpz_set_ui(A->y,0);
	A->Inf=TRUE;
}

void EFp_Clear(struct EFp *A){
	mpz_clear(A->x);
	mpz_clear(A->y);
}

void EFp_Rand(struct EFp *A){
	unsigned int r;
	struct EFp Ans_Copy;
	mpz_t tmp1,tmp2;
	EFp_Init(&Ans_Copy);
	mpz_init(tmp1);
	mpz_init(tmp2);
	gmp_randstate_t state;
	gmp_randinit_default(state);

	while(mpz_legendre(tmp1,prime)!=1){
		r=rand();
		gmp_randseed_ui(state,r);
		mpz_urandomm(Ans_Copy.x,state,prime);
		mpz_powm_ui(tmp1,Ans_Copy.x,3,prime);
		mpz_mul(tmp2,Ans_Copy.x,Ans_Copy.x);
		mpz_mul(tmp2,tmp2,ECCPara_a);
		mpz_add(tmp1,tmp1,tmp2);
		mpz_add(tmp1,tmp1,Ans_Copy.x);
	}
	Fp_Sqrt(Ans_Copy.y,tmp1);
	EFp_Copy(A,&Ans_Copy);

	mpz_clear(tmp1);
	EFp_Clear(&Ans_Copy);
}


void EFp_Show(struct EFp *A){
	gmp_printf("(%Zd,%Zd)\n",A->x,A->y);
}

void EFp_ECA(struct EFp *Ans, struct EFp *P, struct EFp *Q){

	if(P->Inf==TRUE){
		EFp_Copy(Ans, Q);
		return;
	}else if(Q->Inf==TRUE){
		EFp_Copy(Ans, P);
		return;
	}else if(mpz_cmp(P->x,Q->x)==0&&mpz_cmp(P->y,Q->y)){
		EFp_SetInf(Ans);
		return;
	}else if(mpz_cmp(P->x,Q->x)==0&&mpz_cmp(P->y,Q->y)==0){
		EFp_ECD(Ans,P);
		return;
	}

	//lambda = (y1-y2)/(x1-x2)
	//x3 = lambda^2-x1-x2-a
	//y3 = lambda(x1-x3)-y1

	mpz_t lambda,tmp;
	struct EFp Ans_Copy;

	mpz_init(lambda);
	mpz_init(tmp);
	EFp_Init(&Ans_Copy);

	mpz_sub(tmp,P->x,Q->x);
	Fp_Inv(tmp,tmp);
	mpz_sub(lambda,P->y,Q->y);
	mpz_mul(lambda,lambda,tmp);
	mpz_mod(lambda,lambda,prime);

	mpz_mul(Ans_Copy.x,lambda,lambda);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,P->x);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,Q->x);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,ECCPara_a);
	mpz_mod(Ans_Copy.x,Ans_Copy.x,prime);

	mpz_sub(tmp,P->x,Ans_Copy.x);
	mpz_mul(Ans_Copy.y,lambda,tmp);
	mpz_sub(Ans_Copy.y,Ans_Copy.y,P->y);
	mpz_mod(Ans_Copy.y,Ans_Copy.y,prime);

	EFp_Copy(Ans,&Ans_Copy);

	mpz_clear(lambda);
	mpz_clear(tmp);
	EFp_Clear(&Ans_Copy);
}

void EFp_ECD(struct EFp *Ans, struct EFp *P){
	if(P->Inf==TRUE){
		EFp_SetInf(Ans);
		return;
	}else if(mpz_cmp_ui(P->y,0)==0){
		EFp_SetInf(Ans);
		return;
	}

	//lambda = (3x1^2+2ax1+1)/2y1
	//x3 = lambda^2-x1-x2-a
	//y3 = lambda(x1-x3)-y1

	mpz_t lambda,tmp;
	struct EFp Ans_Copy;

	mpz_init(lambda);
	mpz_init(tmp);
	EFp_Init(&Ans_Copy);

	mpz_mul(lambda,P->x,P->x);
	mpz_mul_ui(lambda,lambda,3);
	mpz_mul(tmp,P->x,ECCPara_a);
	mpz_mul_ui(tmp,tmp,2);
	mpz_add(lambda,lambda,tmp);
	mpz_add_ui(lambda,lambda,1);

	mpz_mul_ui(tmp,P->y,2);
	Fp_Inv(tmp,tmp);
	mpz_mul(lambda,lambda,tmp);
	mpz_mod(lambda,lambda,prime);

	mpz_mul(Ans_Copy.x,lambda,lambda);
	mpz_mul_ui(tmp,P->x,2);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,tmp);
	mpz_sub(Ans_Copy.x,Ans_Copy.x,ECCPara_a);
	mpz_mod(Ans_Copy.x,Ans_Copy.x,prime);

	mpz_sub(tmp,P->x,Ans_Copy.x);
	mpz_mul(Ans_Copy.y,lambda,tmp);
	mpz_sub(Ans_Copy.y,Ans_Copy.y,P->y);
	mpz_mod(Ans_Copy.y,Ans_Copy.y,prime);

	EFp_Copy(Ans,&Ans_Copy);
	mpz_clear(lambda);
	mpz_clear(tmp);
	EFp_Clear(&Ans_Copy);
}

void EFp_XZ_Init(struct EFp_XZ *A){
	mpz_init(A->x);
	mpz_init(A->z);
	A->Inf=FALSE;
}

void EFp_XZ_Clear(struct EFp_XZ *A){
	mpz_clear(A->x);
	mpz_clear(A->z);
}

void EFp_XZ_ECD(struct EFp_XZ *Ans, struct EFp_XZ *P){
	mpz_t tmp1,tmp2,tmp3;
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(tmp3);
	struct EFp_XZ Ans_Copy;
	EFp_XZ_Init(&Ans_Copy);

	mpz_mul(tmp1,x,x);
	mpz_mod(tmp1,tmp1,prime);
	mpz_mul(tmp2,z,z);
	mpz_mod(tmp2,tmp2,prime);
	mpz_mul(tmp3,x,z);
	mpz_mod(tmp3,tmp3,prime);
	mpz_mul(Ans_Copy->z,tmp3,ECCPara_a);
	mpz_mod(Ans_Copy->z,Ans_Copy->z,prime);

	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
	EFp_Clear(&Ans_Copy);
}

void EFp_SCM(struct EFp *Ans, struct EFp *P, mpz_t j){
	int i;
	int r;
	r = (int)mpz_sizeinbase(j,2);
	struct EFp Ans_Copy;
	EFp_Init(&Ans_Copy);
	EFp_Copy(&Ans_Copy,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp_ECD(&Ans_Copy,&Ans_Copy);
			EFp_ECA(&Ans_Copy,&Ans_Copy,P);
		}else{
			EFp_ECD(&Ans_Copy,&Ans_Copy);
		}
	}

	EFp_Copy(Ans,&Ans_Copy);
	EFp_Clear(&Ans_Copy);
}

void EFp_XZ_Init(struct EFp_XZ *A){
	
}

//----------------------------------------------------------------
void Set_Parameter(mpz_t prime){
	mpz_init(ECCPara_a);
	mpz_ui_pow_ui(prime,2,255);
	mpz_sub_ui(prime,prime,19);
	mpz_set_ui(ECCPara_a,486662);
}

void hoge(){
	struct EFp A,B,C;
	EFp_Init(&A);
	EFp_Init(&B);
	EFp_Init(&C);

	EFp_Rand(&A);
	EFp_Rand(&B);
	EFp_ECA(&C,&A,&B);
	EFp_ECD(&C,&C);
	EFp_Show(&C);
	EFp_ECD(&A,&A);
	EFp_ECD(&B,&B);
	EFp_ECA(&C,&A,&B);
	EFp_Show(&C);

	EFp_Clear(&A);
	EFp_Clear(&B);
	EFp_Clear(&C);
}

int main(void){

	unsigned int now = (unsigned int)time(0);
	srand(now);
	Set_Parameter(prime);
	hoge();

	return 0;
}
