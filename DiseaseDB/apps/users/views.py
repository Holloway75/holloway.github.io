from django.contrib.auth import authenticate, login, logout
from django.shortcuts import render, redirect
from .forms import CaptchaLoginForm

def login_view(request):
    if request.method == 'POST':
        form = CaptchaLoginForm(data=request.POST)
        if form.is_valid():
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password')
            user = authenticate(request, username=username, password=password)
            if user is not None:
                login(request, user)
                return redirect('home')  # 重定向到主页或其他页面
    else:
        form = CaptchaLoginForm()
    return render(request, 'users/login.html', {'form': form})

def logout_view(request):
    logout(request)
    return redirect('login')  # 登出后重定向到登录页面