# v2,v3 (STAR BES)

Данный код позволяет измерить эллиптический и треугольный потоки в столкновениях ядер золота при энергиях BES STAR, используя метод плоскости события, SP.

## Содержание:

I. [Устновка](#Устновка) \
II. [подготовка данных](#ПодготовкаДанных) \
III. [Обработка событий](#EventProcessing) \
III. [Анализ результатов](#Result) \

## I. Устновка <a name="Устновка"></a>
Требования:
CMake (ver. 3.1 or higher)
ROOT (should work with versions 5 and 6)
ROOT's MathMore library

```bash
mkdir STAR
cd STAR
git clone https://github.com/DemanovAE/FlowAnalysis_new.git
cd FlowAnalysis_new
```

Не забудьте добавить библиотеки ROOT в свою среду, используя thisroot.sh
В терминале:

```sh
source /opt/fairsoft/bmn/may18p1/bin/thisroot.sh
#or
source /mnt/pool/rhic/4/parfenovpeter/Soft/Basov/ROOT/build/bin/thisroot.sh
```

Компилирование читалки данных:
(собранная читалка для данных femtoDst лежит по пути `/scratch2/demanov/STAR/BES/StFemtoEvent/libStFemtoDst.so`)

```bash
cd StFemtoEvent
make
```
Измените пути в `set_env.sh`: изменить `ST_FEMTO_DST_INC_DIR` на стандартный. Это каталог, в котором хранится `libStFemtoDst.so`.

Установка проекта с помощью CMake. (находясь в директории ./macro)

```bash
cd ../
mkdir build
cd build
cmake ../macro
make
```

## II. Подготовка данных <a name="ПодготовкаДанных"></a>

Используйте `./scripts/GenerateLists.sh` для создания списков файлов:

```bash
. GenerateLists.sh FEMTODST_DIR N_FILES_IN_LIST
```
где FEMTODST_DIR - путь к каталогу с femtoDst.root файлами. И N_FILES_IN_LIST обозначает максимальное количество femtoDst.root файлов в каждом листе. Базовый пример: `. GenerateLists.sh /scratch2/parfenov/StData/27gev/run1/ 100`

Полученные списки файлов будут в BES/lists/

В директории уже созданны листы с данными для кластера МИФИ.


## III. Обработка событий <a name="EventProcessing"></a>

Интерактивный режим:
```bash
./FemtoDstAnalyzer_PID -i inFile -o outFile -r inRec -f inFlat -m WorkMode -g Energy
```
1. `inFile` - root файл или лист с файлами
2. `outFile` - выходной файл. Для измерения потоков нужно провести 3 прогонки данных. Выходной файл после 1 прогонки `NoRe_27GeV.root`, после второй - `Re_27GeV.root` и конечный файл `Flow_27GeV.root`.
3. `inRec` -
4. `inFlat` -
5. `WorkMode` - указывает стадию анализа.

        | WorkMode        | Описание    |
        | --------------- | ----------- |
                                                                                                                            |
        | raw             | Первая прогонка данных. На данном этапе набираются данные для корекции Q-векторов, а именно дял процедуры реценренинга              |
        | rec             | Вторая прогонка данных. На данном этапе набираются данные для корекции угла плоскости события, а именно дял процедуры флатенинга    |
        | flow            | Третья прогонка данных. На данном этапе измеряются v2, v3,... и расрешение.                                                            |

6. `Energy` - энергия анализируемых данных.

Пример запуска:
```bash
./FemtoDstAnalyzer_PID -i ../lists/lists27GeV/StRun15.list -o ./raw_27GeV.root -r null -f null -m raw -g 27
./FemtoDstAnalyzer_PID -i ../lists/lists27GeV/StRun15.list -o ./rec_27GeV.root -r raw_27GeV.root -f null -m raw -g 27
./FemtoDstAnalyzer_PID -i ../lists/lists27GeV/StRun15.list -o ./flow_27GeV.root -r raw_27GeV.root -f rec_27GeV.root -m raw -g 27
```

Для отправки задач на класстер используются следующий bash скрипт (!! надо настроить пути в свою директорию) - `/scripts/cl_mephi/start_pid_nica.sh INPUT_FILELIST_DIR INPUT_WORKMODE INPUT_ENERGY`
```sh
cd ../scripts/cl_mephi
. start_pid_nica.sh /scratch2/$USER/STAR/BES/lists/lists27GeV raw 27
```

После каждого прогона данных необходимо объединить все выходные файлы при помощи hadd
```bash
hadd -j 4 ../OUT/27GeV/${WorkMode}_${Energy}GeV.root ../OUT/27GeV/${Energy}_${WorkMode}_*/root/StRun*.root
```

## IIII. Анализ результатов <a name="Result"></a>

Объеденить RunId и удлаить BadRun
```bash
cd /DrawScripts
root -l -b -q FirstScript.cpp"(\"inFileAfterHadd.root\",\"./OutFlowFile.root\","28")"

```

Перейти в директорию DrawScripts

```bash
cd DrawScripts
mkdir build
cd build
cmake ../
make
cd ../PictCode
root -l -b -q PictDiffParAntVsCent.cpp"(\"OutFlowFile.root\","Energy")
```