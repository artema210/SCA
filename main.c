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

struct Fp ECCPara_a;//y^2 = x^3+ax^2+x

struct Fp{
	mpz_t x0;
};

struct EFp{
	struct Fp x,y;
	int Inf;
};

struct EFp_XZ{
	struct Fp x,z;
	int Inf;
};

//----------Fp functions--------------------------------------------
void Fp_Init(struct Fp *A);
void Fp_Set(struct Fp *A, struct Fp *B);
void Fp_Clear(struct Fp *A);
void Fp_Add(struct Fp *Ans, struct Fp *A, struct Fp *B);
void Fp_Mul(struct Fp *Ans, struct Fp *A, struct Fp *B);
void Fp_Pow(struct Fp *Ans, struct Fp *A, struct Fp *B);	//Ans <- A^B
void Fp_Inv(struct Fp *Ans, struct Fp *A);			//Ans <- A^(-1)
void Fp_Div(struct Fp *Ans, struct Fp *A, struct Fp *B);//Ans <- A/B
void Fp_Sqrt(struct Fp *Ans, struct Fp *A);			//ANs <- A^(1/2)
//----------EFp functions--------------------------------------------
void EFp_Init(struct EFp *A);	// vector initialization
void EFp_Set(struct EFp *A, struct EFp *B);	//A <- B
void EFp_SetInf(struct EFp *A);	//A <- Infinity
void EFp_Clear(struct EFp *A);
void EFp_Rand(struct EFp *A);
void EFp_Show(struct EFp *A);
void EFp_Inv(struct EFp *A, struct EFp *B);
void EFp_ECA(struct EFp *Ans, struct EFp *P, struct EFp *Q);//ANS=P+Q
void EFp_ECD(struct EFp *Ans, struct EFp *P);//ANS=2*P
void EFp_SCM(struct EFp *Ans, struct EFp *P,mpz_t j);//ANS=[j]P

void EFp_XZ_Init(struct EFp_XZ *A);
void EFp_XZ_Set(struct EFp_XZ *A, struct EFp_XZ *B);	//A <- B
void EFp_XZ_SetInf(struct EFp_XZ *A);	//A <- Infinity
void EFp_XZ_Clear(struct EFp_XZ *A);
void EFp_Show(struct EFp *A);
void EFp_XZ_ECA(struct EFp_XZ *Ans, struct EFp_XZ *P, struct EFp_XZ *Q, struct Fp *Px);//ECA by XZ 
void EFp_XZ_ECD(struct EFp_XZ *Ans, struct EFp_XZ *P);//ECD by XZ
void EFp_XZ_ML(struct EFp_XZ *Ans, struct EFp_XZ *P, mpz_t j);//SCM by MontgomeryLadder. Ans=[j]P

void EFp_Conv_XZ(struct EFp_XZ *A, struct EFp *B);//convert XY to XZ
void EFp_Conv_XY(struct EFp *A, struct EFp_XZ *B);//convert XZ to XY

//-------------Fp functions--------------------------------
void Fp_Init(struct Fp *A){
	mpz_init(A->x0);
}

void Fp_Set(struct Fp *A, struct Fp *B){
	mpz_set(A->x0, B->x0);
}

void Fp_Clear(struct Fp *A){
	mpz_clear(A->x0);
}

void Fp_Add(struct Fp *Ans, struct Fp *A, struct Fp *B){
	mpz_add(Ans->x0,A->x0,B->x0);
	mpz_mod(Ans->x0,Ans->x0,prime);
}

void Fp_Sub(struct Fp *Ans, struct Fp *A, struct Fp *B){
	mpz_sub(Ans->x0,A->x0,B->x0);
	mpz_mod(Ans->x0,Ans->x0,prime);
}

void Fp_Mul(struct Fp *Ans, struct Fp *A, struct Fp *B){
	mpz_mul(Ans->x0,A->x0,B->x0);
	mpz_mod(Ans->x0,Ans->x0,prime);
}

void Fp_Inv(struct Fp *Ans, struct Fp *A){
	mpz_t j;
	mpz_init(j);
	struct Fp Ans_Copy;
	Fp_Init(&Ans_Copy);
	mpz_init(j);
	Fp_Set(&Ans_Copy,A);
	mpz_sub_ui(j,prime,2);
	int i,r;
	r=(int)mpz_sizeinbase(j,2);
	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			Fp_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);
			Fp_Mul(&Ans_Copy,&Ans_Copy,A);
		}else{
			Fp_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);
		}
	}
	Fp_Set(Ans,&Ans_Copy);

	Fp_Clear(&Ans_Copy);
	mpz_clear(j);
}

void Fp_Pow(struct Fp *Ans, struct Fp *A, struct Fp *B){
	if(mpz_cmp_ui(B->x0,0)==0){
		mpz_set_ui(Ans->x0,1);
		return;
	}
	struct Fp Ans_Copy;
	int i;
	int r;//bit

	Fp_Init(&Ans_Copy);
	r= (int)mpz_sizeinbase(B->x0,2);
	Fp_Set(&Ans_Copy,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B->x0,i)==1){
			Fp_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);	
			Fp_Mul(&Ans_Copy,&Ans_Copy,A);			//(A*2)*A
		}else{
			Fp_Mul(&Ans_Copy,&Ans_Copy,&Ans_Copy);	//A*2
		}
	}
	Fp_Set(Ans,&Ans_Copy);
	Fp_Clear(&Ans_Copy);
}

void Fp_Div(struct Fp *Ans, struct Fp *A, struct Fp *B){
	Fp_Inv(Ans,B);
	Fp_Mul(Ans,Ans,A);
}

void Fp_Sqrt(struct Fp *Ans, struct Fp *A){
	if(mpz_legendre(A->x0,prime)!=1){
		printf("No Answer\n");
		return;
	}

	struct Fp b,e,n,q,r,t,x,y,tmp1,tmp2,base;
	unsigned int m,Rand;
	gmp_randstate_t state;
	gmp_randinit_default(state);

	Fp_Init(&b);
	Fp_Init(&e);
	Fp_Init(&n);
	Fp_Init(&q);
	Fp_Init(&r);
	Fp_Init(&t);
	Fp_Init(&x);
	Fp_Init(&y);
	Fp_Init(&tmp1);
	Fp_Init(&tmp2);
	Fp_Init(&base);
	
	mpz_set_ui(base.x0,2);
	Fp_Set(&n,A);
	while(mpz_legendre(n.x0,prime)!=-1){
		Rand=rand();
		gmp_randseed_ui(state,Rand);
		mpz_urandomm(n.x0, state, prime);
	}
	mpz_set_ui(n.x0,21);
	mpz_set_ui(e.x0,0);
	mpz_sub_ui(q.x0,prime,1);
	
	while(mpz_odd_p(q.x0)==0){
		mpz_add_ui(e.x0,e.x0,1);
		mpz_div_ui(q.x0,q.x0,2);
	}

	Fp_Pow(&y,&n,&q);
	Fp_Set(&r,&e);
	mpz_sub_ui(q.x0,q.x0,1);
	mpz_div_ui(q.x0,q.x0,2);
	Fp_Pow(&x,A,&q);
	Fp_Set(&tmp1,&x);
	Fp_Mul(&x,&x,A);
	Fp_Mul(&b,&x,&tmp1);

	while(mpz_cmp_ui(b.x0,1)){
		m=-1;
		Fp_Init(&tmp2);
		while(mpz_cmp_ui(tmp2.x0,1)){
			m++;
			mpz_ui_pow_ui(tmp1.x0,2,m);
			Fp_Pow(&tmp2,&b,&tmp1);
		}
		mpz_sub_ui(tmp1.x0,r.x0,m);
		mpz_sub_ui(tmp1.x0,tmp1.x0,1);
		Fp_Pow(&tmp2,&base,&tmp1);
		Fp_Pow(&t,&y,&tmp2);
		mpz_powm_ui(y.x0,t.x0,2,prime);
		mpz_set_ui(r.x0,m);
		Fp_Mul(&x,&x,&t);
		Fp_Mul(&b,&b,&y);
	}
	
	Fp_Set(Ans,&x);
	
	Fp_Clear(&b);
	Fp_Clear(&e);
	Fp_Clear(&n);
	Fp_Clear(&q);
	Fp_Clear(&r);
	Fp_Clear(&t);
	Fp_Clear(&x);
	Fp_Clear(&y);
	Fp_Clear(&tmp2);
	Fp_Clear(&tmp1);
	Fp_Clear(&base);

}

//-------------EFp functions--------------------------------
void EFp_Init(struct EFp *A){
	Fp_Init(&A->x);
	Fp_Init(&A->y);
	A->Inf=FALSE;
}

void EFp_Set(struct EFp *A, struct EFp *B){
	Fp_Set(&A->x,&B->x);
	Fp_Set(&A->y,&B->y);
	A->Inf=B->Inf;
}

void EFp_SetInf(struct EFp *A){
	mpz_set_ui(A->x.x0,0);
	mpz_set_ui(A->y.x0,0);
	A->Inf=TRUE;
}

void EFp_Clear(struct EFp *A){
	Fp_Clear(&A->x);
	Fp_Clear(&A->y);
}

void EFp_Rand(struct EFp *A){
	unsigned int r;
	struct EFp Ans_Copy;
	struct Fp tmp1,tmp2;
	EFp_Init(&Ans_Copy);
	Fp_Init(&tmp1);
	Fp_Init(&tmp2);
	gmp_randstate_t state;
	gmp_randinit_default(state);

	while(mpz_legendre(tmp1.x0,prime)!=1){
		r=rand();
		gmp_randseed_ui(state,r);
		mpz_urandomm(Ans_Copy.x.x0,state,prime);
		mpz_powm_ui(tmp1.x0,Ans_Copy.x.x0,3,prime);
		Fp_Mul(&tmp2,&Ans_Copy.x,&Ans_Copy.x);
		Fp_Mul(&tmp2,&tmp2,&ECCPara_a);
		Fp_Add(&tmp1,&tmp1,&tmp2);
		Fp_Add(&tmp1,&tmp1,&Ans_Copy.x);
	}
	Fp_Sqrt(&Ans_Copy.y,&tmp1);
	EFp_Set(A,&Ans_Copy);

	Fp_Clear(&tmp1);
	Fp_Clear(&tmp2);
	EFp_Clear(&Ans_Copy);
}


void EFp_Show(struct EFp *A){
	gmp_printf("(%Zd,%Zd)\n",A->x.x0,A->y.x0);
}

void EFp_Inv(struct EFp *A, struct EFp *B){
	Fp_Set(&A->x,&B->x);
	mpz_neg(A->y.x0,B->y.x0);
	mpz_mod(A->y.x0,A->y.x0,prime);
}

void EFp_ECA(struct EFp *Ans, struct EFp *P, struct EFp *Q){

	if(P->Inf==TRUE){
		EFp_Set(Ans, Q);
		return;
	}else if(Q->Inf==TRUE){
		EFp_Set(Ans, P);
		return;
	}else if(mpz_cmp(P->x.x0,Q->x.x0)==0&&mpz_cmp(P->y.x0,Q->y.x0)){
		EFp_SetInf(Ans);
		return;
	}else if(mpz_cmp(P->x.x0,Q->x.x0)==0&&mpz_cmp(P->y.x0,Q->y.x0)==0){
		EFp_ECD(Ans,P);
		return;
	}

	//lambda = (y1-y2)/(x1-x2)
	//x3 = lambda^2-x1-x2-a
	//y3 = lambda(x1-x3)-y1

	struct Fp lambda,tmp;
	struct EFp Ans_Copy;

	Fp_Init(&lambda);
	Fp_Init(&tmp);
	EFp_Init(&Ans_Copy);

	Fp_Sub(&tmp,&P->x,&Q->x);
	Fp_Inv(&tmp,&tmp);
	Fp_Sub(&lambda,&P->y,&Q->y);
	Fp_Mul(&lambda,&lambda,&tmp);

	Fp_Mul(&Ans_Copy.x,&lambda,&lambda);
	Fp_Sub(&Ans_Copy.x,&Ans_Copy.x,&P->x);
	Fp_Sub(&Ans_Copy.x,&Ans_Copy.x,&Q->x);
	Fp_Sub(&Ans_Copy.x,&Ans_Copy.x,&ECCPara_a);

	Fp_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);

	EFp_Set(Ans,&Ans_Copy);

	Fp_Clear(&lambda);
	Fp_Clear(&tmp);
	EFp_Clear(&Ans_Copy);
}

void EFp_ECD(struct EFp *Ans, struct EFp *P){
	if(P->Inf==TRUE){
		EFp_SetInf(Ans);
		return;
	}else if(mpz_cmp_ui(P->y.x0,0)==0){
		EFp_SetInf(Ans);
		return;
	}

	//lambda = (3x1^2+2ax1+1)/2y1
	//x3 = lambda^2-x1-x2-a
	//y3 = lambda(x1-x3)-y1

	struct Fp lambda,tmp;
	struct EFp Ans_Copy;

	Fp_Init(&lambda);
	Fp_Init(&tmp);
	EFp_Init(&Ans_Copy);

	Fp_Mul(&lambda,&P->x,&P->x);
	mpz_mul_ui(lambda.x0,lambda.x0,3);
	Fp_Mul(&tmp,&P->x,&ECCPara_a);
	mpz_mul_ui(tmp.x0,tmp.x0,2);
	Fp_Add(&lambda,&lambda,&tmp);
	mpz_add_ui(lambda.x0,lambda.x0,1);

	mpz_mul_ui(tmp.x0,P->y.x0,2);
	Fp_Inv(&tmp,&tmp);
	Fp_Mul(&lambda,&lambda,&tmp);

	Fp_Mul(&Ans_Copy.x,&lambda,&lambda);
	mpz_mul_ui(tmp.x0,P->x.x0,2);
	Fp_Sub(&Ans_Copy.x,&Ans_Copy.x,&tmp);
	Fp_Sub(&Ans_Copy.x,&Ans_Copy.x,&ECCPara_a);

	Fp_Sub(&tmp,&P->x,&Ans_Copy.x);
	Fp_Mul(&Ans_Copy.y,&lambda,&tmp);
	Fp_Sub(&Ans_Copy.y,&Ans_Copy.y,&P->y);

	EFp_Set(Ans,&Ans_Copy);
	Fp_Clear(&lambda);
	Fp_Clear(&tmp);
	EFp_Clear(&Ans_Copy);
}

void EFp_SCM(struct EFp *Ans, struct EFp *P, mpz_t j){
	int i;
	int r;
	r = (int)mpz_sizeinbase(j,2);
	struct EFp Ans_Copy;
	EFp_Init(&Ans_Copy);
	EFp_Set(&Ans_Copy,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp_ECD(&Ans_Copy,&Ans_Copy);
			EFp_ECA(&Ans_Copy,&Ans_Copy,P);
		}else{
			EFp_ECD(&Ans_Copy,&Ans_Copy);
		}
	}

	EFp_Set(Ans,&Ans_Copy);
	EFp_Clear(&Ans_Copy);
}

void EFp_XZ_Init(struct EFp_XZ *A){
	mpz_init(A->x.x0);
	mpz_init(A->z.x0);
	A->Inf=FALSE;
}

void EFp_XZ_Set(struct EFp_XZ *A, struct EFp_XZ *B){
	Fp_Set(&A->x,&B->x);
	Fp_Set(&A->z,&B->z);
}

void EFp_XZ_SetInf(struct EFp_XZ *A){
	mpz_set_ui(A->x.x0,1);
	mpz_set_ui(A->z.x0,0);
	A->Inf=TRUE;
}

void EFp_XZ_Clear(struct EFp_XZ *A){
	mpz_clear(A->x.x0);
	mpz_clear(A->z.x0);
}

void EFp_XZ_Show(struct EFp_XZ *A){
	gmp_printf("(%Zd,%Zd)\n",A->x.x0,A->z.x0);
}

void EFp_XZ_ECA(struct EFp_XZ *Ans, struct EFp_XZ *P, struct EFp_XZ *Q, struct Fp *Px){
	/*if(P->Inf==TRUE){
		EFp_XZ_Set(Ans, Q);
		return;
	}else if(Q->Inf==TRUE){
		EFp_XZ_Set(Ans, P);
		return;
	}else if(mpz_cmp(P->x.x0,Q->x.x0)==0&&mpz_cmp(P->z.x0,Q->z.x0)){
		EFp_XZ_SetInf(Ans);
		return;
	}else if(mpz_cmp(P->x.x0,Q->x.x0)==0&&mpz_cmp(P->z.x0,Q->z.x0)==0){
		EFp_XZ_ECD(Ans,P);
		return;
	}*/

	struct Fp tmp;
	Fp_Init(&tmp);
	struct EFp_XZ Ans_Copy;
	EFp_XZ_Init(&Ans_Copy);

	Fp_Mul(&Ans_Copy.x,&P->x,&Q->x);
	Fp_Mul(&tmp,&P->z,&Q->z);
	Fp_Sub(&Ans_Copy.x,&Ans_Copy.x,&tmp);
	Fp_Mul(&Ans_Copy.x,&Ans_Copy.x,&Ans_Copy.x);

	Fp_Mul(&Ans_Copy.z,&P->x,&Q->z);
	Fp_Mul(&tmp,&P->z,&Q->x);
	Fp_Sub(&Ans_Copy.z,&Ans_Copy.z,&tmp);
	Fp_Mul(&Ans_Copy.z,&Ans_Copy.z,&Ans_Copy.z);
	Fp_Mul(&Ans_Copy.z,&Ans_Copy.z,Px);
	EFp_XZ_Set(Ans,&Ans_Copy);

	Fp_Clear(&tmp);
	EFp_XZ_Clear(&Ans_Copy);
}

void EFp_XZ_ECD(struct EFp_XZ *Ans, struct EFp_XZ *P){
	struct Fp tmp1,tmp2,tmp3;
	Fp_Init(&tmp1);
	Fp_Init(&tmp2);
	Fp_Init(&tmp3);
	struct EFp_XZ Ans_Copy;
	EFp_XZ_Init(&Ans_Copy);

	Fp_Mul(&tmp1,&P->x,&P->x);
	Fp_Mul(&tmp2,&P->z,&P->z);
	Fp_Mul(&tmp3,&P->x,&P->z);

	Fp_Sub(&Ans_Copy.x,&tmp1,&tmp2);
	Fp_Mul(&Ans_Copy.x,&Ans_Copy.x,&Ans_Copy.x);

	Fp_Mul(&Ans_Copy.z,&tmp3,&ECCPara_a);
	Fp_Add(&Ans_Copy.z,&Ans_Copy.z,&tmp1);
	Fp_Add(&Ans_Copy.z,&Ans_Copy.z,&tmp2);
	Fp_Mul(&Ans_Copy.z,&Ans_Copy.z,&tmp3);
	mpz_mul_ui(Ans_Copy.z.x0,Ans_Copy.z.x0,4);
	mpz_mod(Ans_Copy.z.x0,Ans_Copy.z.x0,prime);
	
	EFp_XZ_Set(Ans,&Ans_Copy);

	Fp_Clear(&tmp1);
	Fp_Clear(&tmp2);
	Fp_Clear(&tmp3);
	EFp_XZ_Clear(&Ans_Copy);
}

void EFp_XZ_ML(struct EFp_XZ *Ans, struct EFp_XZ *P, mpz_t j){
	int i;
	int r;//num of bits
	r = (int)mpz_sizeinbase(j,2);
	struct EFp_XZ Ans_Copy,tmp;
	EFp_XZ_Init(&Ans_Copy);
	EFp_XZ_Init(&tmp);
	EFp_XZ_Set(&tmp,P);
	EFp_XZ_SetInf(&Ans_Copy);

	for(i=r-1;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp_XZ_ECA(&Ans_Copy,&Ans_Copy,&tmp,&P->x);
			EFp_XZ_ECD(&tmp,&tmp);
		}else{
			EFp_XZ_ECA(&tmp,&Ans_Copy,&tmp,&P->x);
			EFp_XZ_ECD(&Ans_Copy,&Ans_Copy);
		}
	}
	EFp_XZ_Set(Ans,&Ans_Copy);
	EFp_XZ_Clear(&Ans_Copy);
	EFp_XZ_Clear(&tmp);
}

void EFp_Conv_XZ(struct EFp_XZ *A, struct EFp *B){
	Fp_Set(&A->x,&B->x);
	mpz_set_ui(A->z.x0,1);
}

void EFp_Conv_XY(struct EFp *A, struct EFp_XZ *B){
	struct Fp tmp1,tmp2;
	Fp_Init(&tmp1);
	Fp_Init(&tmp2);
	struct EFp Ans_Copy;
	EFp_Init(&Ans_Copy);
	Fp_Div(&Ans_Copy.x,&B->x,&B->z);

	mpz_pow_ui(tmp1.x0,Ans_Copy.x.x0,3);
	mpz_mod(tmp1.x0,tmp1.x0,prime);
	Fp_Mul(&tmp2,&Ans_Copy.x,&Ans_Copy.x);
	Fp_Mul(&tmp2,&tmp2,&ECCPara_a);
	Fp_Add(&tmp1,&tmp1,&tmp2);
	Fp_Add(&tmp1,&tmp1,&Ans_Copy.x);
	Fp_Sqrt(&Ans_Copy.y,&tmp1);
	
	EFp_Set(A,&Ans_Copy);
	
	Fp_Clear(&tmp1);
	Fp_Clear(&tmp2);
	EFp_Clear(&Ans_Copy);
}

//----------------------------------------------------------------
void Set_Parameter(){
	Fp_Init(&ECCPara_a);
	mpz_ui_pow_ui(prime,2,255);
	mpz_sub_ui(prime,prime,19);
	mpz_set_ui(ECCPara_a.x0,486662);
}

void hoge(){
	struct EFp A,B;
	struct EFp_XZ C,D;
	EFp_Init(&A);
	EFp_Init(&B);
	EFp_XZ_Init(&C);
	EFp_XZ_Init(&D);
	mpz_t j;
	mpz_init(j);
	mpz_set_ui(j,13423);

	EFp_Rand(&A);
	EFp_SCM(&B,&A,j);
	EFp_Show(&B);
	EFp_Conv_XZ(&C,&A);
	EFp_XZ_ML(&D,&C,j);

	EFp_Conv_XY(&B,&D);
	EFp_Show(&B);
	EFp_Inv(&B,&B);
	EFp_Show(&B);
}

void hoge2(){
	struct EFp A;
	EFp_Init(&A);
	EFp_Rand(&A);
	EFp_Show(&A);
	EFp_SCM(&A,&A,prime);
	EFp_Show(&A);

}

int main(void){

	unsigned int now = (unsigned int)time(0);
	srand(now);
	Set_Parameter(prime);
	hoge();

	return 0;
}
