const csrftoken = getCookie('csrftoken');
const sendToServer = (id, value, type) => {
    $.ajax({
        url: '/updateOrgdetail/', //"{% url 'organism_updatedetail' %}", 
        type: "POST",
        headers: { 'X-CSRFToken': csrftoken },
        data: { id: id, value: value, type: type }

    })
        .done((response) => {
            console.log(response)
            // window.location.reload()
        })
        .fail(() => {
            console.log("Error occured")
        })
}
