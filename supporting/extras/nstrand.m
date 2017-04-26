def nstrand(u, p = None, return_ef = False):
    '''
    %   Calculate the coeficients of multi-channel auto-regressive matrix using
    %   Nuttall-Strand algorithm (a generalization of single channel harmonic
    %                             method)
    %
    %   Input parameters:
    %     IP     - Ordem of autoregressive model (integer)
    %     u      - Complex matrix with NUMCHS channels of sample data
    %
    %   Output parameters:
    %     PF     - Covariance matrix of NUMCHS x NUMCHS of linear forward
    %              prediction error
    %     A      - Complex array of forward linear prediction matrix
    %              coefficients
    %     PB     - Complex backward linear prediction error covariance array
    %     B      - Complex array of backward linear prediction matrix
    %              coefficients
    '''
    [lx,cx]=u.shape
    if lx > cx:
        print ('Input matrix is probably transposed.')
        return
    
    if p is None:
        p = pr_.maxp
    
    NUMCHS=lx      #% Number of channels.
    MAXORDER=200   #% Maximum order of AR model allowed for calculation.
    N=max(u.shape) #% N - Number of samples per channel.
    IP = p

    #    Initialization
    ISTAT=0
    if (IP > MAXORDER):
        ISTAT=3
        print('IP > 200')
        return
    
    ef=u.copy()                    #% Eq. (15.91)
    eb=u.copy()                    #% Eq. (15.91)
    pf=dot(u, u.transpose())       #% Eq. (15.90)
    pb=array(pf)                   #% Eq. (15.90)
    M=0
    #    Main Loop
    while 1:
        #%  Update estimated covariance errors  Eq. (15.89)
        pfhat=dot(ef[:,M+1:N],ef[:,M+1:N].transpose())
        pbhat=dot(eb[:,M:N-1],eb[:,M:N-1].transpose())
        pfbhat=dot(ef[:,M+1:N],eb[:,M:N-1].transpose())
        M=M+1
        #%  Calculate estimated partial correlation matrix - Eq. (15.98)
        #%             (Nuttall-Strand algorithm only)
        RHO=lyap(dot(pfhat,inv(pf)),dot(inv(pb),pbhat),-2*pfbhat)
        #%  Update forward and backward reflection coeficients
        #%  Eqs. (15.73),(15.74),(15.78) (algoritmo de Nuttall-Strand)
        AM=dot(-RHO,inv(pb))
        BM=dot(-RHO.transpose(),inv(pf))
        dimA=AM.shape[0]
        dimB=BM.shape[0]
        if M == 1:
            A=zeros((1,dimA,dimA),float)
            B=zeros((1,dimB,dimB),float)
        else:
            #print 'M', M
            #print 'A', A
            #print 'AM', AM
            #print 'AA', A.resize((M,dimA,dimA))
            #print 'dimA', dimA
            A=resize(A,(M,dimA,dimA))
            B=resize(B,(M,dimB,dimB))
        A[M-1,:,:]=AM
        B[M-1,:,:]=BM
        #%  Update forward and backward covariance error  - Eqs. (15.75),(15.76)
        pf=pf-dot(dot(AM,BM),pf)
        pb=pb-dot(dot(BM,AM),pb)
        
        #%  Update forward and backward predictor coeficients - Eqs.(15.84),(15.85)
        if not (M == 1):
            for K in range(1,M):
                temp1=A[K-1,:,:].copy()
                A[K-1,:,:]=A[K-1,:,:]+dot(AM,B[M-K-1,:,:])
                B[M-K-1,:,:]=B[M-K-1,:,:]+dot(BM,temp1)
        #%  Update residues
        #%  existe erro no calculo dos residuos
        Tef=array(ef)
        ef[0:NUMCHS,range(N-1,M-1,-1)]=ef[:,range(N-1,M-1,-1)]+dot(AM,eb[:,range(N-2,M-2,-1)])
        eb[0:NUMCHS,range(N-1,M-1,-1)]=eb[:,range(N-2,M-2,-1)]+dot(BM,Tef[:,range(N-1,M-1,-1)])
        
        #%  Verify if model order is adequate
        if M == IP:
            A=-A
            B=-B
            break
    
    #print 'pf:', pf
    #pf = abs(pf) #TODO: conferir o pf, pq as vezes da matriz toda negativa?
    
    if (return_ef):
        return A.transpose(1,2,0), ef #TODO: porque o abs?
    else:
        #return pf,A,pb,B,ef,eb,ISTAT 
        return A.transpose(1,2,0), pf/N