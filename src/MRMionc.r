library(parallel)


args = commandArgs(trailingOnly=TRUE)

if (length(args)>0) {
    wstr=""
} else{
    wstr="2"
}
trans_files=Sys.glob(file.path('trans',paste0("trans",wstr,"_*Q1_*Q3_*.txt")))
dir.create('IonChromatogram')

mzML_files=(readLines("trans_mzML_list.txt"))

colorp=c('blue','green','red','cyan','magenta')
colorp=rep(colorp,2)

by_sample=0

if (by_sample ){
    dat=vector(length=length(trans_files))
    for (i in 1:length(trans_files)) {
        dat[i]=list(readLines(trans_files[i]))
        print(trans_files[i])
    }
    pos=rep(2,length(trans_files))

    for (i in seq(1,length(mzML_files),10000)) {
        pdf(file.path('IonChromatogram',paste0('IonC',wstr,'_',mzML_files[i],'.pdf')))#,useDingBats = FALSE)
        par(mfrow=c(4,2),oma = c(0,0,0,0) , mar = c(2,2,2,1))
        for (ii in 1:length(trans_files)) {
            dat0=unlist(dat[ii])
            rt_=unlist(strsplit(dat0[pos[ii]],"\t"))
            I_=as.numeric(unlist(strsplit(dat0[pos[ii]+1],"\t")))
            peaks=c()
            for (iii in (pos[ii]+2):length(dat0)) {
                if (dat0[iii]!=""){
                    peaks=c(peaks,list(as.numeric(unlist(strsplit(dat0[iii],"\t")))))
                }else{
                    pos[ii]=iii+1
                    break
                }
            }
            plot(x=rt_,y=I_,pch=19,cex=.3,xlab='',ylab='',ylim=c(0,max(I_)*1.1))
            title(main = dat0[1])
            minI=min(I_)
            if(length(peaks)>0)for (iii in 1:min(length(peaks),length(colorp))){
                peak=unlist(peaks[iii])
                if (wstr==''){
                    abline(v=peak[1],col=colorp[iii])
                    points(x=c(peak[3],peak[4]),y=rep(minI,2),col=colorp[iii])
                } else {
                    abline(v=peak[1],col=colorp[iii],lty=peak[3])
                    peaklen=(length(peak)-3)/2
                    points(x=peak[3+1:peaklen],y=peak[peaklen+3+1:peaklen],col=colorp[iii],cex=.7)
                }
            }
        }
        dev.off()
        if (i%%100==0) print(mzML_files[i])
    }

}

fx <- function(trans_file) {
    dat=readLines(trans_file)
    if (length(dat)<=2){
        print('empty')
        quit()
    }
    startdat=c(1)
    for (i in 1:length(dat)){
        if (dat[i]==''){
            startdat=c(startdat,i)
        }
    }
    find1_=unlist(gregexpr('_',trans_file))[1]
    filename_=paste0('IonC',substr(trans_file,find1_,nchar(trans_file)-4),'.pdf')
    pdf(file.path('IonChromatogram',filename_))#,useDingBats = FALSE)
    print(filename_)
    par(mfrow=c(4,2),oma = c(0,0,0,0) , mar = c(2,2,2,1))
    for (i in 1:(length(startdat)-1)){#for each transition
        start_=startdat[i]+1
        end_=startdat[i+1]-1
        rt_=as.numeric(unlist(strsplit(dat[start_],"\t")))
        I_=as.numeric(unlist(strsplit(dat[start_+1],"\t")))
        peaks=c()
        if (start_+1<end_)for (ii in (start_+2):end_){
            peaks=c(peaks,list(as.numeric(unlist(strsplit(dat[ii],"\t")))))
        }

        if (rt_[length(rt_)]-rt_[1]>200 && length(peaks)>0) {
            x_min=unlist(peaks[1])[1]-60
            x_max=unlist(peaks[length(peaks)])[1]+60
            plot(x=rt_,y=I_,pch=19,cex=.3,xlab='',ylab='',ylim=c(0,max(I_)*1.1),xlim=c(x_min,x_max))
        } else {
            plot(x=rt_,y=I_,pch=19,cex=.3,xlab='',ylab='',ylim=c(0,max(I_)*1.1))
        }

        title(main = substr(mzML_files[i],0,nchar(mzML_files[i])-5))#, sub = "Sub-title",
        minI=min(I_)
        if(length(peaks)>0)for (i in 1:min(length(peaks),length(colorp))){
            peak=unlist(peaks[i])
            if (wstr==''){
                abline(v=peak[1],col=colorp[i])
                points(x=c(peak[3],peak[4]),y=rep(minI,2),col=colorp[i])
            } else {
                abline(v=peak[1],col=colorp[i],lty=peak[3])
                peaklen=(length(peak)-3)/2
                points(x=peak[3+1:peaklen],y=peak[peaklen+3+1:peaklen],col=colorp[i],cex=.7)
            }
        }
    }
    dev.off()
}

results = mclapply(trans_files, fx, mc.cores =19)


fx <- function(trans_file) {
    dat=readLines(trans_file)
    if (length(dat)<=2){
        print('empty')
        quit()
    }
    startdat=c(1)
    for (i in 1:length(dat)){
        if (dat[i]==''){
            startdat=c(startdat,i)
        }
    }
    find1_=unlist(gregexpr('_',trans_file))[1]

    max_I_=0
    for (i in 1:(length(startdat)-1)){#for each transition
        if (grepl('TQC',mzML_files[i])) {
            start_=startdat[i]+1
            end_=startdat[i+1]-1
            I_=as.numeric(unlist(strsplit(dat[start_+1],"\t")))
            max_I_=max(max_I_,I_)
        }
    }

    filename_=paste0('TQC',substr(trans_file,find1_,nchar(trans_file)-4),'.png')
    png(file.path('IonChromatogram',filename_))
    firstflag=0
    for (i in 1:(length(startdat)-1)){#for each transition
        if (grepl('TQC',mzML_files[i])) {
            start_=startdat[i]+1
            end_=startdat[i+1]-1
            rt_=as.numeric(unlist(strsplit(dat[start_],"\t")))
            I_=as.numeric(unlist(strsplit(dat[start_+1],"\t")))
            if (firstflag==0) {
                plot(x=rt_,y=I_,xlab='RT',ylab='intensity',ylim=c(0,max_I_)*1.0,type='l',col=rgb(red=0.0, green=0.0, blue=0.0, alpha=0.2))
                firstflag=firstflag+1
            } else {
                lines(x=rt_,y=I_,col=rgb(red=0.0, green=0.0, blue=0.0, alpha=0.2))
                firstflag=firstflag+1
            }
        }
    }
    print(paste(filename_,firstflag))
    dev.off()
}


results = mclapply(trans_files, fx, mc.cores =19)

