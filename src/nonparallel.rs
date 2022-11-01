use std::{f32::consts::PI, ops::Index};
use image::{DynamicImage, GenericImageView};
use ndarray::parallel::prelude::{IntoParallelIterator, ParallelIterator};
use num::{traits::real, Float};

const THREADS: usize = 4;


fn double_syn(x: i32, y: i32) -> f32{
    let x_f = x as f32;
    let y_f = y as f32;
    return x_f.sin() + y_f.sin();
}


fn generate_numbers(x: Vec<i32>, y: Vec<i32>, f: &dyn Fn(i32, i32) -> f32) -> Vec<Vec<f32>>{
    let mut res: Vec<Vec<f32>> = vec!();
    for x_ in x{
        let mut cur_res: Vec<f32> = vec!();
        for y_ in &y{
            cur_res.push(f(x_, y_.to_owned()));
        }
        res.push(cur_res);
    }
    return res;
}

fn create_test_image(size: i32, f: &dyn Fn(i32, i32) -> f32) -> ndarray::Array2::<num::complex::Complex32>{
    let x: Vec<i32> = (0..size).collect();
    let y: Vec<i32> = (0..size).collect();
    let img = vec_to_ndarray(generate_numbers(x, y, f));
    return img;
}


fn create_matrix(size: usize) -> ndarray::Array2::<num::complex::Complex32>{
    let mut F = ndarray::Array2::<num::complex::Complex32>::default((size, size));
    for (j, mut row) in F.axis_iter_mut(ndarray::Axis(0)).enumerate(){
        for (k, col) in row.iter_mut().enumerate(){
            *col = exp_matrix(j as f32, k as f32, size as f32);
        }
    }
    return F;
}

fn vec_to_ndarray(input: Vec<Vec<f32>>) -> ndarray::Array2::<num::complex::Complex32>{
    let shape = (input.len(), input[0].len());
    let mut res = ndarray::Array2::<num::complex::Complex32>::default(shape);
    for (i, mut row) in res.axis_iter_mut(ndarray::Axis(0)).enumerate() {
        for (j, col) in row.iter_mut().enumerate() {
            *col = num::complex::Complex32::new(input[i][j], 0.0);
        }
    }
    return res;
}

fn exp_matrix(j: f32, k: f32, N: f32) -> num::complex::Complex32{
    let mut real_part = 2.0 * PI * j * k / N;
    real_part = real_part.cos();
    let mut im_part = -2.0 * PI * j * k / N;
    im_part = im_part.sin();
    let res = num::complex::Complex::new(real_part, im_part);
    return res;
}

fn dft(input: &ndarray::Array2::<num::complex::Complex32>) -> ndarray::Array2::<num::complex::Complex32>{
    let shape = input.shape();
    let M = shape[0];
    let N = shape[1];
    let F_M = create_matrix(M as usize);
    let F_N = create_matrix(N as usize);
    let mut res = F_M.dot(input);
    res = res.dot(&F_N);
    return res;
}


fn image_to_vec(img: DynamicImage) -> Vec<Vec<f32>>{
    let mut res: Vec<Vec<f32>> = vec!();
    let dimensions = img.dimensions();
    for i in 0..(dimensions.0 as usize){
        res.push(vec!());
        for j in 0..dimensions.1{
            res[i].push(img.get_pixel(i as u32, j).0[0] as f32);
        } 
    }
    return res;
}

fn prep_magnitude(input: &ndarray::Array2::<num::complex::Complex32>) -> (ndarray::Array2::<f32>, f32, f32){
    let shape = (input.shape()[0], input.shape()[1]);
    let mut res = ndarray::Array2::<f32>::default(shape);
    let mut max: f32 = 0.0;
    let mut min: f32 = 0.0;
    for (i, mut row) in res.axis_iter_mut(ndarray::Axis(0)).enumerate() {
        for (j, col) in row.iter_mut().enumerate() {
            let curr = (1.0 + input[[i, j]].norm()).log2();
            if curr > max{
                max = curr;
            } else if curr < min{
                min = curr;
            }
            *col = curr;
        }
    }
    return (res, max, min);
}

fn normalize_magnitude(min: f32, max: f32, value: f32) ->f32{
    let res = (value - min)* (255.0/(max - min));
    return res;

}


fn process_image_gen(image: &ndarray::Array2::<num::complex::Complex32>) -> f32{
    let start = std::time::Instant::now();
    let test_dft = dft(&image);
    let duration = start.elapsed();
    return duration.as_nanos() as f32 / 1000000.0;
}

fn save_processed_image(img: &Vec<Vec<f32>>, res_img:  &mut image::GrayImage, input: &ndarray::Array2::<num::complex::Complex32>, path: &str){
    let (mag_dft, max, min )= prep_magnitude(&input);
    for (r, row) in img.iter().enumerate(){
       for (c, col) in row.iter().enumerate(){
           let pix_val = *col as u8;
           let pix_val_mgntd = normalize_magnitude(min, max, mag_dft[[r, c]]) as u8;
           res_img.put_pixel(r.try_into().unwrap(), c.try_into().unwrap(), image::Luma([pix_val]));
           res_img.put_pixel(
            (img.len() + (r + img.len()/2)%(img.len())).try_into().unwrap(),
            ((c + img[0].len()/2)%(img[0].len())).try_into().unwrap(),
           
         image::Luma([pix_val_mgntd]));
       }
    }
    let res_path = "res_".to_owned() + path;
    res_img.save(res_path).unwrap();
}

fn process_image(path: &str){
    let img = image_to_vec(image::open(path).unwrap().grayscale());
    let mut res_img: image::GrayImage = image::ImageBuffer::new((img.len()*2).try_into().unwrap(), img[0].len().try_into().unwrap());
    let img_nd = vec_to_ndarray(img.clone());
    let test_dft = dft(&img_nd);
    save_processed_image(&img, &mut res_img, &test_dft, path);
}

fn sample_dft(size: i32, sample_size: usize, f: &dyn Fn(i32, i32) -> f32){
    let test_img = create_test_image(size, f);
    let mut total_time = 0.0;
    for i in 0..sample_size{
        let curr_time = process_image_gen(&test_img);
        total_time += curr_time;
    }
    println!("Mean processing time {}", total_time / (sample_size as f32));
}


fn main2(){
    //debug and see images
    let path = "230.jpg";
    process_image(path);
}

fn main() {
    let size = std::env::args().nth(1).expect("size must be unsigned int").parse::<i32>().unwrap();
    let sample_size = std::env::args().nth(2).expect("sample size must be unsigned int").parse::<usize>().unwrap();
    sample_dft(size, sample_size, &double_syn);
}
