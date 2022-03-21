# 作者：韩开元
# 单位：内蒙古大学 物理科学与技术学院 2018级 数理基础科学(基地，物理学)
# 版本：内测版
# 最后更新时间：2021 年 12 月 16 日 19:08

# 目标：创建一个 SmileTools 父类，预计在其下创建 bed 子类，smile 子类分别用来分析处理 bed 文件和 smile 文件

# 最大的父类，包含对于 bed 文件或 smile 文件的通用操作函数
class SmileTools:

    # 此父类仅为包含通用函数，一般不作为数据类型
    def __init__(self):
        self.list = []
        self.annotation = "Don\'t build SmileTools data，please use BedFile or SmileFile format data."
        self.type = "SmileTools"

    # 将 a_list 中每一行用 '\t' 分隔后认为每行 type_index 元素为类型，展示这个列表中存在哪些类型
    def which_type(self, type_index):
        types = []
        for line in self.list:
            type_i = line.split('\t')[type_index]
            if type_i not in types:
                types.append(type_i)
        return types

    # 将 find_list 列表中每一行按照 '\t' 分隔后，找出 type_index 索引的元素等于 name 的行
    def find_type_i(self, type_index, name):
        find_list = self.list
        result_list = []
        for line in find_list:
            line_already_split = line.split('\t')
            type_i = line_already_split[type_index]
            if type_i == name:
                result_list.append(line)
        return result_list

    # 在 find_list 中找出包含 name 的行
    def query(self, name):
        find_list = self.list
        result = []
        for line in find_list:
            if name in line:
                result.append(line)
        return result

    # 将 a_list 的每一行用 '\t 分隔后，展示第 index 索引的信息
    def show_index(self,  index, a_list=[]):
        if type(self) == BedFile:
            a_list = self.list
        for i in range(len(a_list)):
            show = a_list[i].split('\t')[index]
            print([i, show])

    # 用 '\t' 分隔符分隔每一行然后展示出来
    def show_by_split(self, a_list=[], show_range='all'):
        if type(self) == BedFile:
            a_list = self.list
        if show_range == 'all':
            for line in a_list:
                print(line.split('\t'))
        else:
            for line in a_list[show_range[0] : show_range[1]]:
                print(line.split('\t'))

    # 将 ori_list 用 '\t' 分隔后，按照第 index 索引处的数据大小，将整个列表排序 —— 仅为 sort_index 服务
    def sort_index_for_sort_index(self, ori_list, index):
        demo_list = []
        for line in ori_list:
            demo_list.append(line)

        def find_index(a_line):
            left = int(a_line.split('\t')[index])
            return left

        demo_list.sort(key=find_index)
        return demo_list

    # 对 bed/smile 的数据，对于每个 type 内的数据按照 sort_index 索引处数据的大小进行排序
    def sort_index(self, type_index, sort_index):
        result_loops = []
        types = self.which_type(type_index)
        for type_i in types:
            temp_loops = self.find_type_i(type_index, type_i)
            mid_loops = self.sort_index_for_sort_index(temp_loops, sort_index)
            for loop in mid_loops:
                result_loops.append(loop)
        return result_loops


# BedFile 子类，是对所有 bed 类型文件的处理
class BedFile(SmileTools):

    # 包含：列表类型、列表注释、列表本体 三个成员变量
    def __init__(self, origin_list=[], list_type='Unknown'):
        super().__init__()
        self.type = list_type
        annotation_list = []
        result_list = []
        for line in origin_list:
            if len(line) > 3:
                if line[0] != '#':
                    if line[0:3] == 'chr':
                        result_list.append(line[3:])
                    else:
                        result_list.append(line)
                else:
                    annotation_list.append(line)
        self.annotations = annotation_list
        self.list = result_list

    # 从指定路径读取文件，返回 BedFile 文件
    def bed_load_path(path, list_type):
        f = open(path, 'r')
        read_list = f.read().split('\n')
        f.close()
        return BedFile(read_list, list_type)

    # 将 BedFile 文件按照 chr_index 索引为染色体，在染色体内按 sort_index 位置进行排序
    def chr_sort(self, chr_index, sort_index):
        result_loops = self.sort_index(type_index=chr_index, sort_index=sort_index)
        return BedFile(origin_list=result_loops, list_type=self.type)

    # 根据两个元素 cluster_index 之间距离阈值 distance_threshold 来聚类，输出一个新的 BedFile
    # 在最后增加四列：类序列、类距离平均值、与上一个元素的距离、与下一个元素的距离
    def cluster_by_distance(self, distance_threshold, cluster_indexes=[1, 2]):
        # 将返回一个新的 BedFile 文件
        new_bed_file = BedFile(origin_list=self.list, list_type=self.type)
        # 提取出 start 和 end 索引
        start_index = cluster_indexes[0]
        end_index = cluster_indexes[1]
        # 每一个元素属于第几个簇，第一个元素一定是第一个簇内的
        in_same = [1]
        distance_with_last = []
        distance_with_next = []
        average_distance = []
        # 第 j 个簇
        j = 1
        # # 判断当前元素与上一个元素是否是在同一簇内，0 为不是，1 为是，认为第一个元素是在簇内的！
        if_in = 1
        # 用 for 循环完成聚类
        for i in range(len(new_bed_file.list) - 1):
            # 识别当前与下一个元素的染色体与始末位置
            chromosome_1 = new_bed_file.list[i].split('\t')[0]
            chromosome_2 = new_bed_file.list[i+1].split('\t')[0]
            # 始末位置用 int 转换成整数方便计算距离
            # 找到上一个元素的终止位置
            if i == 0:
                last_end = int(new_bed_file.list[i].split('\t')[start_index])
            else:
                last_end = int(new_bed_file.list[i-1].split('\t')[end_index])
            # 当前元素的起止位置
            start = int(new_bed_file.list[i].split('\t')[start_index])
            end = int(new_bed_file.list[i].split('\t')[end_index])
            # 下一个元素的起始位置
            next_start = int(new_bed_file.list[i+1].split('\t')[start_index])
            # 计算距离并添加至列表中
            distance_last = start - last_end
            distance_next = next_start - end
            distance_with_last.append(distance_last)
            distance_with_next.append(distance_next)
            # 如果与下一个元素的距离在阈值 distance_threshold 内：
            if chromosome_1 == chromosome_2 and distance_next <= distance_threshold:
                # 如果开始了一个新的簇，簇类加一
                if if_in == 0:
                    if_in = 1
                # 如果还在同一个簇内没有新的操作
            # 与下一个元素不在同一个阈值 distance_threshold 内：
            else:
                # 结束了当前簇 且 簇类加一
                if_in = 0
                j = j + 1
            # 将下一个元素所属簇加入 in_same 列表
            in_same.append(j)
        # 将最后一个元素起止增添进去
        temp_last_end = int(new_bed_file.list[len(new_bed_file.list)-2].split('\t')[end_index])
        temp_start = int(new_bed_file.list[len(new_bed_file.list)-1].split('\t')[start_index])
        temp_distance_last = temp_start - temp_last_end
        distance_with_last.append(temp_distance_last)
        distance_with_next.append(0)
        # 用来记录每个簇内元素间距
        cluster_distance = []
        # 中间值使用的列表
        mid_distance = []
        # 计算均值
        for i in range(len(in_same)-1):
            # 在一个簇内时，计算距离并存在 mid_distance 内
            if in_same[i + 1] == in_same[i]:
                mid_distance.append(distance_with_next[i])
            # 不在一个簇内时，将 mid_distance 增添进 cluster_distance 内并归零 mid_distance
            else:
                cluster_distance.append(mid_distance)
                mid_distance = []
        # 最后一次循环之后最后的 mid_distance 没有添加进去，需要手动添加！
        cluster_distance.append(mid_distance)

        # 定义求均值函数
        def average_int(a_list):
            sum_i = 0
            for num in a_list:
                sum_i = sum_i + num
            return int(sum_i/len(a_list))

        # 计算每个 簇 的均值，只有一个元素的簇返回 '.'
        for i in in_same:
            # 簇类从 1 开始，索引从 0 开始，注意要 -1 ！
            if len(cluster_distance[i-1]) > 0:
                average_distance.append(average_int(cluster_distance[i-1]))
            else:
                average_distance.append('.')
        # 将上述元素增添至新 BedFile 文件
        # 最后增添四个元素：簇类、簇内平均距离、与上个元素距离、与下个元素距离
        for i in range(len(new_bed_file.list)):
            new_line = new_bed_file.list[i] + "\t" + str(in_same[i]) + "\t" + str(average_distance[i]) + "\t"
            new_line = new_line + str(distance_with_last[i]) + "\t" + str(distance_with_next[i])
            new_bed_file.list[i] = new_line
        return new_bed_file

    # 将一个 BedFile 全部类型数据存储至 out_path 路径处的名为 file_name 的文件
    def save_as(self, out_folder, file_name="handled.bed", out_indexes="all"):
        out_path = out_folder + file_name
        fff = open(out_path, 'w+')
        fff.close()
        fw = open(out_path, 'a')
        if out_indexes == "all":
            for i in range(len(self.list)):
                data = self.list[i] + '\n'
                fw.write(data)
        else:
            for i in range(len(self.list)):
                mid_line = self.list[i].split('\t')
                new_line = ''
                for j in range(len(out_indexes) - 1):
                    new_line = new_line + mid_line[out_indexes[j]] + '\t'
                new_line = new_line + mid_line[out_indexes[len(out_indexes)-1]] + '\n'
                fw.write(new_line)
        fw.close()
        print("Complete～～～")

    # 搜寻某个 name 在 chromosome 染色体上的范围：
    # 若是 gene_name 可以在 name 处用严格搜索即："\"gene_name\""
    def query_range(self, name, location_indexes=[1, 2]):
        query_list = self.query(name=name)
        query_bed_file = BedFile(origin_list=query_list, list_type=self.type)
        chromosomes = query_bed_file.which_type(type_index=0)
        range_result = []
        list_result = []
        for chromosome in chromosomes:
            find_list = query_bed_file.find_type_i(type_index=0, name=chromosome)
            list_result.append(find_list)
            left = int(find_list[0].split('\t')[location_indexes[0]])
            right = int(find_list[0].split('\t')[location_indexes[1]])
            for line in find_list:
                left_i = int(line.split('\t')[location_indexes[0]])
                right_i = int(line.split('\t')[location_indexes[1]])
                if left_i < left:
                    left = left_i
                if right_i > right:
                    right = right_i
            mid_result = [chromosome, [left, right]]
            range_result.append(mid_result)
        if len(range_result) == 1:
            return [range_result[0], list_result[0]]
        return [range_result, list_result]

    # 寻找某个范围附近的元素，用 location_index 索引找到该元素的起止位置
    # query_range 应该是 [chromosome, [start_location, end_location]] 形式
    def query_neighbours(self, query_range, limit_distance, location_indexes):
        result = []
        target_chromosome = query_range[0]
        target_start = int(query_range[1][0]) - int(limit_distance)
        target_end = int(query_range[1][1]) + int(limit_distance)
        for ori_line in self.list:
            line = ori_line.split('\t')
            chromosome = line[0]
            start = int(line[location_indexes[0]])
            end = int(line[location_indexes[1]])
            if chromosome == target_chromosome:
                if target_start <= start <= target_end or target_start <= end <= target_end:
                    result.append(ori_line)
                elif start <= target_start <= end or start <= target_end <= end:
                    result.append(ori_line)
        return result

    # 寻找距离某个 query_range 最近的元素，搜寻步长为 step_size
    def query_nearst_project(self, query_range, step_size, location_indexes):
        test = self.query_neighbours(query_range, 0, location_indexes)
        if len(test) > 0:
            print("The nearest project have returned and the nearst distance is 0 ~~~")
            return test
        nearest = 0
        step_i = step_size
        loop_time = 1
        while True:
            result = self.query_neighbours(query_range, nearest, location_indexes)
            # 加速搜寻速度
            if len(result) == 1:
                start = int(result[0].split('\t')[location_indexes[0]])
                end = int(result[0].split('\t')[location_indexes[1]])
                range_start = query_range[1][0]
                range_end = query_range[1][1]
                temp_a = abs(start - range_start)
                temp_b = abs(start - range_end)
                temp_c = abs(end - range_start)
                temp_d = abs(end - range_end)
                nearest_distance = min(temp_a, temp_b, temp_c, temp_d)
                print("The nearest project have returned and the nearst distance is "+str(nearest_distance)+" ~~~")
                return result
            # 防止过度搜索且停止加速，使用二分法最速找到答案
            elif len(result) > 1:
                if step_i > 0:
                    nearest = nearest - step_i
                    step_i = int(step_i / 2)
                    loop_time = 0
                else:
                    print("Sorry, but this our best result ~~~")
                    return result
            # 反向二分法加速搜索速度
            if loop_time == 1:
                step_i = step_i * 2
            nearest = nearest + step_i

    # 输入 gene_name，基因 BedFile，看看该基因是否在 loop 上或搜索阈值内 loops
    # 请尽量使用 query_range + query_neighbours + query_nearest_project ，而少使用此功能
    def query_gene_neighbour_loops(self, gene_name, gene_bed_file, limit_distance, location_indexes):
        if self.type != 'loop':
            print("Your self data have to be \'loop\' type BedFile!")
            return -1
        query_range = gene_bed_file.query_range(gene_name, [1, 2])[0]
        if type(query_range[0]) != str:
            print("Pleas make sure your gene only on a single chromosome!")
            return -1
        result = self.query_neighbours(query_range, limit_distance, location_indexes)
        if len(result) == 0:
            print("There are no loop in your limit distance!")
            print("So we try to find nearest loop to your gene and return the result ~~~")
            result = self.query_nearst_project(query_range, 100, location_indexes)
        return result
        
    # 转换不同版本的基因坐标
    # 使用此功能需要先
    def hg_to_hg(self, origin_hg, to_hg, location_indexes, chromosome_index=0):
        # 导入坐标转换包
        from pyliftover import LiftOver
        lo = LiftOver(origin_hg, to_hg)
        ori_list = self.list
        result_list = []
        for ori_line in ori_list:
            line = ori_line.split('\t')
            ori_chromosome = line[chromosome_index]
            # LiftOver 只能使用 chr 类型的染色体，在此进行检查转换
            if ori_chromosome[0:3] != 'chr':
                chromosome = 'chr' + ori_chromosome
            else:
                chromosome = ori_chromosome
            new_locations = []
            # 对于每个需要转换的坐标进行转换
            for index in location_indexes:
                # LiftOver 需要输入 int 类型的坐标
                location = int(line[index])
                new_location = lo.convert_coordinate(chromosome, location)
                # 判断是否转换成功，可能存在新 hg 中没有当前坐标情况
                if type(new_location) == list and new_location != []:
                    # LiftOver 会输出 int 类型的坐标，需要转换成 str 再添加进我们的数据
                    new_locations.append(str(new_location[0][1]))
            # 判断是否此行的所有坐标都转换成功，若失败则直接跳到下一行
            if len(new_locations) != len(location_indexes):
                continue
            # 若成功，将此行转换为新行
            for i in range(len(location_indexes)):
                line[location_indexes[i]] = new_locations[i]
            new_line = ''
            for i in range(len(line) - 1):
                new_line = new_line + line[i] + '\t'
            new_line = new_line + line[len(line) - 1]
            result_list.append(new_line)
        # 返回一个新的 BedFile 文件
        new_bed_file = BedFile(origin_list=result_list, list_type=self.type)
        return new_bed_file


# SmileFile 子类，是对所有 smile 文件的处理，且只能从 loops 文件得到 SmileFile 文件！
class SmileFile(SmileTools):

    def __init__(self, loops_list=["^-^"]):
        super().__init__()
        self.loop = loops_list
        loops = []
        left_anchor = []
        right_anchor = []
        middle = []
        for i in range(len(loops_list)):
            left_anchor.append([])
            right_anchor.append([])
            middle.append([])
        self.left = left_anchor
        self.right = right_anchor
        self.middle = middle
        self.contain = ['loop']
        self.annotation = ["This is a SmileFile formate data，it can contain chromatin loop and other annotations."]

    # 将一个 loop 类型的 BedFile 读取成初始 smile 文件
    def smile_load_loop(bed_loops):
        result_list = []
        if bed_loops.type == 'loop':
            for line in bed_loops.list:
                result_list.append("^-^\tloop\t" + line)
        return SmileFile(loops_list=result_list)

    def smile_load_TAD(bed_loops):
        result_list = []
        if bed_loops.type == 'TAD':
            for line in bed_loops.list:
                result_list.append("^-^\tTAD\t" + line)
        return SmileFile(loops_list=result_list)

    # 从 smile 文件直接读取至一个 SmileFile
    def smile_load_path(path):
        f = open(path, 'r')
        read_list = f.read().split('\n')
        f.close()
        loops = []
        left = []
        right = []
        middle = []
        annotation = []
        temp_bed_file = BedFile(origin_list=read_list)
        contain = temp_bed_file.which_type(type_index=1)
        loop_index = -1
        for line in read_list:
            if len(line) > 3:
                if line[0] == '<':
                    left[loop_index].append(line)
                elif line[0] == '=':
                    middle[loop_index].append(line)
                elif line[0] == '>':
                    right[loop_index].append(line)
                elif line[0:3] == '^-^':
                    loops.append(line)
                    loop_index = loop_index + 1
                    left.append([])
                    right.append([])
                    middle.append([])
                elif line[0] == '#':
                    annotation.append(line)
        bed_loops = BedFile(origin_list=loops, list_type='loop')
        smiles = SmileFile()
        smiles.loop = loops
        smiles.left = left
        smiles.right = right
        smiles.middle = middle
        smiles.annotation = annotation
        smiles.contain = contain
        return smiles

    # 将一个 BedFile 文件增添至一个 SmileFile
    # add 函数中 indexes 函数表示该 BedFile 代表位置的索引，建议将初始文件更改为默认格式
    def add(self, bed_file, indexes=[1,2], anchor_plus=0):
        # 将 smile 文件中的 loop 增添至一个 BedFile 类型的 loops 列表
        loops = BedFile(self.loop, 'loop')
        # 默认已经拥有当下类型注释
        if_contained = 1
        # 判断此 smile 包含类型中是否已经有当下类型了，若没有，则添加
        if bed_file.type not in self.contain:
            self.contain.append(bed_file.type)
            # 第一次添加该类型数据，if_contained 变为 0
            if_contained = 0
        # 找出 loops 存在哪些染色体
        chromosomes = loops.which_type(type_index=2)
        # 对于每一个染色体
        for chromosome in chromosomes:
            # 找到 smile 文件中处于该染色体上所有 loop
            chr_i_loops = loops.find_type_i(type_index=2, name=chromosome)
            # 找到 bed 文件中处于该染色体上所有 annotation
            chr_i_annotation = bed_file.find_type_i(type_index=0, name=chromosome)
            # 对于此染色体上所有 loop 做遍历，将每个 loop 上的 annotation 找出
            for ori_loop in chr_i_loops:
                # 找到这个 loop 在整个 smile.loop 中的索引
                loop_num = self.loop.index(ori_loop)
                # 记得用 '\t' 分隔每一行
                loop = ori_loop.split('\t')
                # 将这个 loop 左右 anchor 的范围分别找出，根据要求增宽 anchor 范围
                left_anchor_range = [int(loop[3])-anchor_plus, int(loop[4])+anchor_plus]
                right_anchor_range = [int(loop[6])-anchor_plus, int(loop[7])+anchor_plus]
                # 对于此染色体上所以 annotation 做遍历，找出处于该 loop 上的 annotation
                for ori_annotation in chr_i_annotation:
                    # 记得用 '\t' 分隔每一行
                    annotation = ori_annotation.split('\t')
                    # 对于所有能代表此 annotation 位置的索引数据
                    for index in indexes:
                        # 用整型转换得到代表该 annotation 位置的 location
                        location = int(annotation[index])
                        # 如果这个 location 在 left_anchor 上
                        if left_anchor_range[0] <= location <= left_anchor_range[1]:
                            # 用 smile 格式处理此行 annotation
                            new_line = ("<\t" + bed_file.type + "\t" + ori_annotation)
                            # 判断此注释是否已经在 smile 文件中了！不可重复输入相同注释！
                            if new_line not in self.left[loop_num]:
                                # 未重复时：找到此 loop 对应的 smile.left 并将 annotation 加入！
                                self.left[loop_num].append(new_line)
                        # 同上，判断是否处于 right_anchor 上
                        elif right_anchor_range[0] <= location <= right_anchor_range[1]:
                            new_line = (">\t" + bed_file.type + "\t" + ori_annotation)
                            if new_line not in self.right[loop_num]:
                                self.right[loop_num].append(new_line)
                        # 同上，判断是否处于 loop 上
                        elif left_anchor_range[1] <= location <= right_anchor_range[0]:
                            new_line = ("=\t" + bed_file.type + "\t" + ori_annotation)
                            if new_line not in self.middle[loop_num]:
                                self.middle[loop_num].append(new_line)
                # 在增加完所有此 loop 拥有的 annotation 以后，整理排序
                # 若这不是第一次增添此类型数据！
                if if_contained == 1:
                    left_not_sorted = BedFile(origin_list=self.left[loop_num])
                    right_not_sorted = BedFile(origin_list=self.right[loop_num])
                    middle_not_sorted = BedFile(origin_list=self.middle[loop_num])
                    self.left[loop_num] = left_not_sorted.sort_index(type_index=1, sort_index=3)
                    self.right[loop_num] = right_not_sorted.sort_index(type_index=1, sort_index=3)
                    self.middle[loop_num] = middle_not_sorted.sort_index(type_index=1, sort_index=3)

    # 将一个 SmileFile 中指定 out_contain 类型存储至 out_path 路径处的名为 file_name 的文件
    def save_as(self, out_folder, file_name="smile.smile", out_contain="all", out_type="anchors"):
        out_path = out_folder + file_name
        fff = open(out_path, 'w+')
        fff.close()
        fw = open(out_path, 'a')
        # 如果是全部写入
        if out_contain == "all":
            for i in range(len(self.loop)):
                data = self.loop[i] + '\n'
                if out_type == "anchors":
                    for line in self.left[i]:
                        data = data + line + '\n'
                    for line in self.right[i]:
                        data = data + line + '\n'
                elif out_type == "all":
                    for line in self.left[i]:
                        data = data + line + '\n'
                    for line in self.middle[i]:
                        data = data + line + '\n'
                    for line in self.right[i]:
                        data = data + line + '\n'
                fw.write(data)
        # 如果是选择性写入...，但是希望永远不要用到这端代码！
        else:
            for i in range(len(self.loop)):
                data = self.loop[i] + '\n'
                if out_type == "anchors":
                    for line in self.left[i]:
                        if line.split('\t')[1] in out_contain:
                            data = data + line + '\n'
                    for line in self.right[i]:
                        if line.split('\t')[1] in out_contain:
                            data = data + line + '\n'
                elif out_type == "all":
                    for line in self.left[i]:
                        if line.split('\t')[1] in out_contain:
                            data = data + line + '\n'
                    for line in self.middle[i]:
                        if line.split('\t')[1] in out_contain:
                            data = data + line + '\n'
                    for line in self.right[i]:
                        if line.split('\t')[1] in out_contain:
                            data = data + line + '\n'
                fw.write(data)
        fw.close()
        print("Complete～～～")

    # 将 a_list 中每一行用 '\t' 分隔后认为每行 type_index 元素为类型，展示这个列表中存在哪些类型
    def which_type(self, type_index):
        types = []
        find_list = self.loop
        for line in find_list:
            type_i = line.split('\t')[type_index]
            if type_i not in types:
                types.append(type_i)
        return types

