// This script is a hacky way to make some small html customizations to the book
$(document).ready(function() {
    // Point the IonSim logo to ionsim.org
    var logolink = document.querySelector('img.logo').parentElement;
    logolink.href = 'https://www.ionsim.org';

    // Remove 'this book' from searchbar placeholder
    $('#search-input').attr('placeholder', 'search examples...');

    // Update copyright year
    var now = new Date();
    var year = now.getFullYear();
    var footertext = document.querySelector('.footer > .container > p');
    if (year < 2029) {
        footertext.innerText = `Copyright © ${year} ionsim.org`;
    }
    // footertext.parentElement.style.backgroundColor = "#003262";
    // footertext.style.color = "white";
    // footertext.style.margin = 0;
    // footertext.style.paddingBottom = "1rem";
    // footertext.style.paddingTop = "1rem";
    // footertext.style.textAlign = "center";
    // var topper = document.querySelector('div.topbar-main');
    // topper.style.backgroundColor = "#003262";
    // var topperbtn = document.querySelector('button.topbarbtn')
    // topperbtn.style.backgroundColor = "transparent"
});