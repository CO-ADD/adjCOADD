$(document).ready(function () {
    $('option').each(function () {
        console.log("selectt")
        replace_string = $(this).text().replace("<small class='not-visible'>", " | ").replace("</small>", "")
        $(this).text(replace_string)
    })
})