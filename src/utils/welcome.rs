pub fn print_welcome() {
    let banner = r#"
    __                           _       
   / /_   _____  ____ _  _   __ (_) ___ 
  / __ \ / ___/ / __ `/ | | / // // _ \
 / /_/ // /    / /_/ /  | |/ // //  __/
/_.___//_/     \__,_/   |___//_/ \___/ 
                                    
Bravie: A (work in progress) Rust DFT Software
----------------------------------------------
    "#;
    println!("{}", banner);
}
