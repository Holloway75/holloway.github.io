from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import CustomUser

class CustomUserAdmin(UserAdmin):
    # 添加自定义字段到管理界面
    list_display = ('username', 'name', 'hospital', 'is_staff')
    fieldsets = UserAdmin.fieldsets + (
        ('自定义信息', {'fields': ('name', 'phone', 'hospital')}),
    )
    add_fieldsets = UserAdmin.add_fieldsets + (
        ('自定义信息', {'fields': ('phone', 'hospital')}),
    )

admin.site.register(CustomUser, CustomUserAdmin)
