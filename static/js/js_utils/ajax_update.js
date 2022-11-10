const csrftoken = getCookie('csrftoken');
const sendToServer = (id, value, type) => {

    console.log(id, typeof (value), type)
    // strain_type_value=

    $.ajax({
        url: '/aa_chem/updateOrgdetail/', //"{% url 'organism_updatedetail' %}", 
        type: "POST",
        headers: { 'X-CSRFToken': csrftoken },
        data: { id: id, value: value, type: type }

    })
        .done((response) => {
            console.log(response)
            window.location.reload()
        })
        .fail(() => {
            console.log("Error occured")
        })
}
