当前脚本用于临时分析空转capture_virus


step1.生成space(fq to bam)、capture_virus(bam to mtx) 分析命令
step2.投递sjm_space.sjm
step3.处理 space bam, 转成 capture_virus 需要的 virus_Aligned.out.bam(query_id 为barcode_umi)
step4.投递sjm_virus.sjm
step5.映射病毒 umi


最优解是基于 celescope 调整分析流程。