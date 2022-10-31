import subprocess
import tqdm
import matplotlib.pyplot as plt

PATH = "target/release/dft-course.exe"

def main():
    sample_size = 10
    sizes = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
    result_times = []
    for size in tqdm.tqdm(sizes):
        command = PATH + " " + str(size) + " " + str(sample_size)
        compl_process = subprocess.run(command.split(' '), capture_output=True)
        temp = compl_process.stdout.decode('utf8').strip()
        result_times.append(float(temp.split(' ')[-1]))
    plt.title("Дискретное преобразование Фурье, Rust")
    plt.xlabel("Размер изображения (size x size)")
    plt.ylabel("Время выполнения, мс")
    plt.scatter(sizes, result_times)
    plt.plot(sizes, result_times)
    plt.show()
main()