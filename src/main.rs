fn generate_numbers(x: Vec<i32>, y: Vec<i32>, f: &dyn Fn(&i32, &i32) -> i32) -> Vec<i32>{
    let mut res: Vec<i32> = vec!();
    for (x_, y_) in x.iter().zip(y){
        res.push(f(&x_, &y_));
    }
    return vec!();
}

fn main() {

}
