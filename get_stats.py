import subprocess
import tqdm
import matplotlib.pyplot as plt

PATH_1 = "target/release/parallel"
PATH_2 = "target/release/nonparallel"

def main():
    sample_size = 10
    sizes = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600]
    result_times_parallel = []
    result_time_nonparallel = []

    print("running parallel...")
    for size in tqdm.tqdm(sizes):
        command1 = PATH_1 + " " + str(size) + " " + str(sample_size)
        compl_process = subprocess.run(command1.split(' '), capture_output=True)
        temp = compl_process.stdout.decode('utf8').strip()
        result_times_parallel.append(float(temp.split(' ')[-1]))

    print("running non-parallel...")
    for size in tqdm.tqdm(sizes):
        command2 = PATH_2 + " " + str(size) + " " + str(sample_size)
        compl_process = subprocess.run(command2.split(' '), capture_output=True)
        temp = compl_process.stdout.decode('utf8').strip()
        result_time_nonparallel.append(float(temp.split(' ')[-1]))

    plt.title("Дискретное преобразование Фурье, Rust")
    plt.xlabel("Размер изображения (size x size)")
    plt.ylabel("Время выполнения, мс")

    plt.scatter(sizes, result_times_parallel)
    plt.plot(sizes, result_times_parallel, label = "Параллельное")

    plt.scatter(sizes, result_time_nonparallel)
    plt.plot(sizes, result_time_nonparallel, label = "Последовательное")
    plt.legend()
    print("parallel ",result_times_parallel)
    print("nonparallel ", result_time_nonparallel)
    plt.show()
main()