use std::f32::consts::PI;
use image::{DynamicImage, GenericImageView};
use rand::Rng;
use num::traits::real;

fn double_syn(x: i32, y: i32) -> f32{
    let x_f = x as f32;
    let y_f = y as f32;
    let mut rng = rand::thread_rng();
    let n1: f32 = rng.gen();
    //return n1;
    return 122.0 * x_f.sin() + 122.0 * y_f.sin();
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

fn exp_euler(k: f32, l: f32, m: f32, n: f32, M: f32, N: f32) -> num::complex::Complex32{
    let mut real_part = 0.0;
    let mut im_part = 0.0;
    real_part = 2.0 * PI * (((k * m / M) + (l * n /N)) as f32);
    real_part = real_part.cos();
    im_part = 2.0 * PI * (((k * m / M) + (l * n /N)) as f32);
    im_part = -1.0 * im_part.sin();
    let res = num::complex::Complex::new(real_part, im_part);
    return res;
}

fn dft_single_value(input: &Vec<Vec<f32>>, k: usize, l: usize, M: f32, N: f32) -> num::complex::Complex32{
    let mut res = num::complex::Complex::new(0.0, 0.0);
    for (m, row) in input.iter().enumerate(){
        let mut cur_res = num::complex::Complex::new(0.0, 0.0);
        for (n, col) in row.iter().enumerate(){
            cur_res += col * exp_euler(
                k as f32,
                l as f32, 
                m as f32, 
                n as f32, 
                M as f32, 
                N as f32);
        }
        res += cur_res;
    }
    return res;
}

fn dft(input: &Vec<Vec<f32>>) -> Vec<Vec<num::complex::Complex32>>{
    let M = input.len() as f32;
    let N = input[0].len() as f32;
    let factorMN = 1.0 / (M * N);
    let mut res: Vec<Vec<num::complex::Complex32>> = vec![vec![num::complex::Complex::new(0.0, 0.0); input.len()]; input[0].len()];
    for (k, row) in input.iter().enumerate(){
        println!("k {}", k);
        for (l, col) in row.iter().enumerate(){
            println!("l {}", l);
            res[k][l] = factorMN * dft_single_value(&input, k, l, M, N);
        }
    }
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

fn main() {
 let x: Vec<i32> = (0..100).collect();
 let y: Vec<i32> = (0..50).collect();
 let mut img_orig: image::GrayImage = image::ImageBuffer::new(x.len().try_into().unwrap(), y.len().try_into().unwrap());
 let mut img_dft_real: image::GrayImage = image::ImageBuffer::new(x.len().try_into().unwrap(), y.len().try_into().unwrap());
 let mut img_dft_im: image::GrayImage = image::ImageBuffer::new(x.len().try_into().unwrap(), y.len().try_into().unwrap());
 let test = generate_numbers(x, y, &double_syn);
 //let img = image_to_vec(image::open("G1-OC4ljwCc.jpg").unwrap().grayscale());
 let test_dft = dft(&test);
 for (r, row) in test.iter().enumerate(){
    for (c, col) in row.iter().enumerate(){
        let pix_val = (*col * 255.0) as u8;
        let pix_val_real = (test_dft[r][c].re * 255.0) as u8;
        let pix_val_im = (test_dft[r][c].im * 255.0) as u8;
        img_orig.put_pixel(r.try_into().unwrap(), c.try_into().unwrap(), image::Luma([pix_val]));
        img_dft_real.put_pixel(r.try_into().unwrap(), c.try_into().unwrap(), image::Luma([pix_val_real]));
        img_dft_im.put_pixel(r.try_into().unwrap(), c.try_into().unwrap(), image::Luma([pix_val_im]));
    }
 }

 img_orig.save("orig.png").unwrap();
 img_dft_real.save("real.png").unwrap();
 img_dft_im.save("im.png").unwrap();
}
