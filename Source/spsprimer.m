function f = spsprimer(Seqlist,Seqname,spectemp,Xratio,Xthreshold,...
    Mgconc,primlength,mismatchparameter,revprim,Cr,Ct)

for i = 1:size(Seqlist,2)
    Seqlistc{i} = seqcomplement(Seqlist{i});
end

SeqArray = [Seqlist;Seqlistc];
mataddress = 0;
Xt = @(G,T,Cr,Ct) 0.5.*((1+((Cr+exp((G.*1000)./...
            (1.987.*(273.15+T))))/Ct))-sqrt((1+((Cr+exp((G.*1000)./...
            (1.987.*(273.15+T))))/Ct)).^2-(4.*(Cr/Ct))));

dGrevprim = ThermoDhyb(revprim,seqrcomplement(revprim),spectemp,Mgconc);

whilestop = 0;
while whilestop == 0
    Candidateprim = {0};
    Candidateprimmat = [];
    msg = sprintf('Selecting %dbp primers for %s against',primlength, Seqname{1});
    xbar=waitbar(0,'Setting up parallel pool');
    titleHandle = get(findobj(xbar,'Type','axes'),'Title');
    set(titleHandle,'FontSize',8)

    
    %Sequence list iterator
    for m = 2:size(SeqArray,2)
        waitbar((m-2)/(size(SeqArray,2)-1),xbar,strcat(msg,...
            {' '}, Seqname{m},' (',...
            string(floor((m-2)/(size(SeqArray,2)-1)*100)),'%)'));

        %Sequence iterator
        parfor w = 2:length(SeqArray{1})-primlength            
            Seq1 = SeqArray{1}(w:w+primlength-1);
            Seq2 = SeqArray{2,m}(w-1:w+primlength);
            numbermismatch = mismatchno(Seq1,seqcomplement(Seq2));
            if numbermismatch >= mismatchparameter                         %    mismatch presence condition
                nonconsenmatseq1 = [Seq1=='A';Seq1=='T';Seq1=='C';Seq1=='G'];
                nonconsenmatseq2 = [Seq2=='A';Seq2=='T';Seq2=='C';Seq2=='G'];
                if sum(sum(nonconsenmatseq1)==0) == 0 && sum(sum(nonconsenmatseq2)==0) == 0         %   intraspecific variable site condition
                    hairpinpres1 = rnafold(Seq1);
                    if sum(find(hairpinpres1 == '(')) == 0                                  
                        %Estimate minimum free energy values and corresponding structures
                        hybridpropm = ThermoDhyb(Seq1,seqrcomplement(SeqArray{1}(w-1:w+primlength)),spectemp,Mgconc);
                        hybridpropmm = ThermoDhyb(Seq1,seqreverse(Seq2),spectemp,Mgconc);
                        Xm = Xt(hybridpropm,spectemp,Cr,Ct);
                        Xmm = Xt(hybridpropmm,spectemp,Cr,Ct);
                        Xmomm = Xm/Xmm;
                        if Xmomm >= Xratio && Xm >= Xthreshold
                            CandidateprimSeq(w,:) = Seq1;
                            CandidateprimHyb(w,:) = hybridpropmm;
                            CandidateprimAdd(w,:) = findstr(SeqArray{1},Seq1);
                        end
                    end
                    
                end

            end
             
        end
        CandidateprimSeq = char(CandidateprimSeq(CandidateprimAdd ~=0,:));
        CandidateprimHyb = CandidateprimHyb(CandidateprimAdd ~=0,:);
        CandidateprimAdd = CandidateprimAdd(CandidateprimAdd ~=0,:);
        
        Candidateprim{m-1,1} = CandidateprimSeq;
        Candidateprim{m-1,2} = CandidateprimHyb;
        Candidateprim{m-1,3} = int64(CandidateprimAdd);
        
        clear CandidateprimSeq CandidateprimHyb CandidateprimAdd

    end
    
    waitbar(1,xbar,strcat(msg,...
            {' '}, Seqname{m},' (100%)'))
    
    if sum(~cellfun(@isempty,Candidateprim)) == length(SeqArray)-1
        Candidateprimmat = union(Candidateprim{1},Candidateprim{2},'rows');
        for i = 1:length(SeqArray)-1
            Candidateprimmat = intersect(Candidateprimmat,Candidateprim{i},'rows');
        end
        cond2 = isempty(Candidateprimmat);
        if cond2 == 0
            whilestop = 1;
            fprintf('\n')
            fprintf('%d',size(Candidateprimmat,1))
            fprintf(' species-specific primer(s) successfully passed the set of criteria')
            fprintf('\n')
            fprintf('\n')
        else
            fprintf('\n')
            fprintf('no primers passed the set of criteria against ')
            fprintf(Seqname{m}(1:end))
            fprintf(',')
            fprintf('\n')
            fprintf('repeating the process with decreased primer length')
            fprintf('\n')
            fprintf('\n')
        end
    end
    
    
    for p = 1:size(Candidateprimmat,1)
        for q = 1:size(Candidateprim,1)
            deltaGarray{p}(q) = Candidateprim{q,2}(strmatch(Candidateprimmat(p,:),Candidateprim{q,1}));
            posarray(p) = Candidateprim{q,3}(strmatch(Candidateprimmat(p,:),Candidateprim{q,1}));
        end
    end
    primlength = primlength-1;
    
    close(xbar)
end

% constructing array for candidate primers and properties
deltaGarray = deltaGarray';
onemat = ones(1,size(Candidateprimmat,1));
Primarray = mat2cell(Candidateprimmat,onemat);

excelarray = {};
excelarray(2:1+length(Primarray),1) = Primarray;
for i = 1:size(deltaGarray,1)
    for j = 2:size(deltaGarray{1},2)+1
        excelarray{i+1,j} = deltaGarray{i}(j-1);
    end
end

for i = 1:size(deltaGarray,1)
    excelarray{i+1,size(deltaGarray{1},2)+2} = posarray(i);
end

for j = 2:size(deltaGarray{1},2)+1
    excelarray{1,j} = Seqname{j}(1:end);                                   %    change range according to sample name in fasta file
end
excelarray{1} = ['Candidate Primer Sequence'];
excelarray{1,end} = 'Primer binding position';

f = excelarray;

end
