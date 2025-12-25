more  PROJECT_above1G_1.txt | awk -F " " '$8 == 2021 || $8==2022 || $8==2023 || $8 ==2024 {print $0}' | grep '\.bam$' > del_bam.txt
more  PROJECT_above1G_1.txt | awk -F " " '$8 == 2021 || $8==2022 || $8==2023 || $8 ==2024 {print $0}' | grep '\.fq$' > del_fq.txt
more  PROJECT_above1G_1.txt | awk -F " " '$8 == 2021 || $8==2022 || $8==2023 || $8 ==2024 {print $0}' | grep '\.tmp$' > del_tmp.txt
more  PROJECT_above1G_1.txt | awk -F " " '$8 == 2021 || $8==2022 || $8==2023 || $8 ==2024 {print $0}' | grep '\.temp$' > del_temp.txt